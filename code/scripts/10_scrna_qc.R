#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(scater)
  library(scDblFinder)
  library(SoupX)
  library(Seurat)
  library(BiocParallel)
  library(jsonlite)
  library(ggplot2)
  library(patchwork)
})

source("scripts/lib_common.R")
theme_set(theme_nature_contract())

args <- commandArgs(trailingOnly = TRUE)
raw_root <- if (length(args) >= 1) args[[1]] else raw_data_root()
processed_root <- if (length(args) >= 2) args[[2]] else "data/processed"
results_root <- if (length(args) >= 3) args[[3]] else "results"
logs_root <- if (length(args) >= 4) args[[4]] else "logs"

dir.create(file.path(processed_root, "scrna", "per_sample"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_root, "tables", "qc"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_root, "figures", "supplementary"), recursive = TRUE, showWarnings = FALSE)
dir.create(logs_root, recursive = TRUE, showWarnings = FALSE)

sample_manifest <- utils::read.delim("metadata/audit/scrna_samples.tsv", check.names = FALSE)
sample_manifest$matrix_dir <- file.path(
  raw_root, "scrna", "matrices",
  paste0(sample_manifest$patient_id, "_", sample_manifest$side, "_raw_feature_bc_matrix")
)

if (!all(dir.exists(sample_manifest$matrix_dir))) {
  missing <- sample_manifest$matrix_dir[!dir.exists(sample_manifest$matrix_dir)]
  stop("Missing matrix directories: ", paste(missing, collapse = ", "))
}

set.seed(722100)
bp <- SerialParam(progressbar = FALSE)
qc_summaries <- list()
threshold_rows <- list()
ambient_rows <- list()
cell_metadata <- list()
seurat_objects <- list()

process_one <- function(row) {
  sample_id <- paste(row$patient_id, row$side, sep = "__")
  message("Reading raw droplets for ", sample_id)
  sce_raw <- DropletUtils::read10xCounts(row$matrix_dir, col.names = TRUE, BPPARAM = bp)
  symbols <- rowData(sce_raw)$Symbol
  symbols[is.na(symbols) | symbols == ""] <- rownames(sce_raw)[is.na(symbols) | symbols == ""]
  symbols <- make.unique(symbols)
  rownames(sce_raw) <- symbols

  message("Calling cells with EmptyDrops for ", sample_id)
  br <- DropletUtils::barcodeRanks(counts(sce_raw), lower = 100)
  ed <- DropletUtils::emptyDrops(counts(sce_raw), lower = 100, niters = 30000,
                                 test.ambient = TRUE, BPPARAM = bp)
  called <- which(!is.na(ed$FDR) & ed$FDR <= 0.01)
  if (length(called) < 200) {
    stop("EmptyDrops called fewer than 200 cells for ", sample_id, "; manual review required")
  }
  sce <- sce_raw[, called]
  colData(sce)$emptydrops_total <- ed$Total[called]
  colData(sce)$emptydrops_logprob <- ed$LogProb[called]
  colData(sce)$emptydrops_fdr <- ed$FDR[called]
  colData(sce)$sample_id <- sample_id
  colData(sce)$patient_id <- row$patient_id
  colData(sce)$side <- row$side
  colData(sce)$age <- row$age
  colData(sce)$sex <- row$sex

  is_mito <- grepl("^MT-", rownames(sce), ignore.case = FALSE)
  is_ribo <- grepl("^RP[SL][0-9]", rownames(sce), ignore.case = FALSE)
  qc <- scater::perCellQCMetrics(sce, subsets = list(Mito = is_mito, Ribo = is_ribo))
  colData(sce) <- cbind(colData(sce), qc)
  thresholds <- robust_qc_thresholds(qc$sum, qc$detected, qc$subsets_Mito_percent)
  passes_metrics <- qc$sum >= thresholds$min_counts & qc$sum <= thresholds$max_counts &
    qc$detected >= thresholds$min_features & qc$detected <= thresholds$max_features &
    qc$subsets_Mito_percent <= thresholds$max_percent_mt
  colData(sce)$passes_qc_metrics <- passes_metrics
  colData(sce)$scDblFinder.score <- NA_real_
  colData(sce)$scDblFinder.class <- "not_assessed_low_quality"

  message("Running scDblFinder for ", sample_id)
  sce_candidate <- scDblFinder::scDblFinder(sce[, passes_metrics], BPPARAM = bp, verbose = FALSE)
  colData(sce)$scDblFinder.score[passes_metrics] <- colData(sce_candidate)$scDblFinder.score
  colData(sce)$scDblFinder.class[passes_metrics] <- as.character(colData(sce_candidate)$scDblFinder.class)
  passes_qc <- passes_metrics & colData(sce)$scDblFinder.class == "singlet"
  colData(sce)$passes_qc <- passes_qc

  message("Estimating ambient RNA with SoupX for ", sample_id)
  toc <- counts(sce)[, passes_qc, drop = FALSE]
  quick <- Seurat::CreateSeuratObject(toc, project = sample_id, min.cells = 0, min.features = 0)
  quick <- Seurat::NormalizeData(quick, verbose = FALSE)
  quick <- Seurat::FindVariableFeatures(quick, nfeatures = min(2000, nrow(quick)), verbose = FALSE)
  quick <- Seurat::ScaleData(quick, features = VariableFeatures(quick), verbose = FALSE)
  quick <- Seurat::RunPCA(quick, features = VariableFeatures(quick), npcs = 20, verbose = FALSE)
  quick <- Seurat::FindNeighbors(quick, dims = 1:15, verbose = FALSE)
  quick <- Seurat::FindClusters(quick, resolution = 0.4, verbose = FALSE)
  clusters <- setNames(as.character(Idents(quick)), colnames(quick))

  soup_status <- "success"
  rho <- rep(NA_real_, ncol(toc))
  names(rho) <- colnames(toc)
  adjusted <- toc
  tryCatch({
    channel <- SoupX::SoupChannel(tod = counts(sce_raw), toc = toc, calcSoupProfile = TRUE)
    channel <- SoupX::setClusters(channel, clusters)
    channel <- SoupX::autoEstCont(channel, doPlot = FALSE, forceAccept = TRUE)
    adjusted <- SoupX::adjustCounts(channel, roundToInt = TRUE)
    rho <- channel$metaData[colnames(toc), "rho"]
  }, error = function(error) {
    soup_status <<- paste("failed", conditionMessage(error), sep = ": ")
  })
  if (soup_status != "success") {
    stop("SoupX correction failed for ", sample_id, ": ", soup_status,
         ". Uncorrected counts must not enter the corrected-count analysis.")
  }

  sce_pass <- sce[, passes_qc]
  assays(sce_pass)$SoupX_counts <- adjusted[, colnames(sce_pass), drop = FALSE]
  colData(sce_pass)$SoupX_rho <- rho[colnames(sce_pass)]
  saveRDS(sce, file.path(processed_root, "scrna", "per_sample", paste0(sample_id, "__called_qc_flags.rds")), compress = FALSE)
  saveRDS(sce_pass, file.path(processed_root, "scrna", "per_sample", paste0(sample_id, "__singlets_soupx.rds")), compress = FALSE)

  meta <- as.data.frame(colData(sce))
  meta$barcode <- rownames(meta)
  meta$cell_id <- paste(sample_id, meta$barcode, sep = "__")
  meta$sample_id <- sample_id

  object <- Seurat::CreateSeuratObject(
    counts = counts(sce_pass), project = sample_id, assay = "RNA",
    meta.data = as.data.frame(colData(sce_pass)), min.cells = 0, min.features = 0
  )
  object[["SoupX"]] <- SeuratObject::CreateAssay5Object(counts = assay(sce_pass, "SoupX_counts"))
  object <- Seurat::RenameCells(object, add.cell.id = sample_id)

  list(
    object = object,
    metadata = meta,
    summary = data.frame(
      sample_id = sample_id, patient_id = row$patient_id, side = row$side,
      raw_barcodes = ncol(sce_raw), emptydrops_called = ncol(sce),
      pass_qc_metrics = sum(passes_metrics), doublets = sum(colData(sce)$scDblFinder.class == "doublet", na.rm = TRUE),
      final_singlets = sum(passes_qc), knee = metadata(br)$knee, inflection = metadata(br)$inflection,
      ambient_status = soup_status, median_soupx_rho = stats::median(rho, na.rm = TRUE)
    ),
    thresholds = data.frame(sample_id = sample_id, patient_id = row$patient_id, side = row$side,
      min_counts = thresholds$min_counts, max_counts = thresholds$max_counts,
      min_features = thresholds$min_features, max_features = thresholds$max_features,
      max_percent_mt = thresholds$max_percent_mt),
    ambient = data.frame(sample_id = sample_id, ambient_status = soup_status,
      median_rho = stats::median(rho, na.rm = TRUE), mean_rho = mean(rho, na.rm = TRUE))
  )
}

for (index in seq_len(nrow(sample_manifest))) {
  result <- process_one(sample_manifest[index, , drop = FALSE])
  sample_id <- result$summary$sample_id
  seurat_objects[[sample_id]] <- result$object
  cell_metadata[[sample_id]] <- result$metadata
  qc_summaries[[sample_id]] <- result$summary
  threshold_rows[[sample_id]] <- result$thresholds
  ambient_rows[[sample_id]] <- result$ambient
  gc(verbose = FALSE)
}

qc_summary <- do.call(rbind, qc_summaries)
threshold_table <- do.call(rbind, threshold_rows)
ambient_table <- do.call(rbind, ambient_rows)
metadata_all <- do.call(rbind, cell_metadata)
write_tsv(qc_summary, file.path(results_root, "tables", "qc", "scrna_qc_summary.tsv"))
write_tsv(threshold_table, file.path(results_root, "tables", "qc", "scrna_qc_thresholds.tsv"))
write_tsv(ambient_table, file.path(results_root, "tables", "qc", "scrna_ambient_summary.tsv"))
write_tsv(metadata_all, file.path(results_root, "tables", "qc", "scrna_cell_qc_metadata.tsv.gz"))

message("Merging QC-passing singlets without expression integration")
combined <- Reduce(function(x, y) merge(x, y), seurat_objects)
DefaultAssay(combined) <- "SoupX"
combined <- JoinLayers(combined, assay = "SoupX")
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, nfeatures = 3000, selection.method = "vst", verbose = FALSE)
combined <- ScaleData(combined, features = VariableFeatures(combined), verbose = FALSE)
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = 30, verbose = FALSE)
combined <- FindNeighbors(combined, dims = 1:25, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:25, reduction = "pca", reduction.name = "umap_unintegrated", verbose = FALSE)
saveRDS(combined, file.path(processed_root, "scrna", "scrna_qc_unintegrated.rds"), compress = FALSE)

umap <- as.data.frame(Embeddings(combined, "umap_unintegrated"))
umap$cell_id <- rownames(umap)
umap$sample_id <- combined$sample_id
umap$patient_id <- combined$patient_id
umap$side <- combined$side
umap$cluster <- as.character(Idents(combined))
write_tsv(umap, file.path(results_root, "tables", "qc", "scrna_unintegrated_umap.tsv"))

plot_meta <- metadata_all
plot_meta$side <- factor(plot_meta$side, levels = c("convex", "concave"))
p1 <- ggplot(plot_meta, aes(x = sample_id, y = log10(sum + 1), fill = side)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_hline(data = threshold_table, aes(yintercept = log10(min_counts + 1)), linetype = 2, linewidth = 0.25) +
  scale_fill_manual(values = palette_contract[c("convex", "concave")]) +
  labs(x = NULL, y = expression(log[10](UMI + 1)), title = "Library complexity") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
p2 <- ggplot(plot_meta, aes(x = sample_id, y = detected, fill = side)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_hline(data = threshold_table, aes(yintercept = min_features), linetype = 2, linewidth = 0.25) +
  scale_fill_manual(values = palette_contract[c("convex", "concave")]) +
  labs(x = NULL, y = "Detected genes", title = "Detected features") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
p3 <- ggplot(plot_meta, aes(x = sample_id, y = subsets_Mito_percent, fill = side)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_hline(data = threshold_table, aes(yintercept = max_percent_mt), linetype = 2, linewidth = 0.25) +
  scale_fill_manual(values = palette_contract[c("convex", "concave")]) +
  labs(x = NULL, y = "Mitochondrial reads (%)", title = "Mitochondrial fraction") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
p4 <- ggplot(qc_summary, aes(x = side, y = final_singlets, group = patient_id, colour = patient_id)) +
  geom_line(linewidth = 0.45) + geom_point(size = 1.5) +
  scale_x_discrete(limits = c("convex", "concave")) +
  labs(x = NULL, y = "QC-passing singlets", title = "Paired retained cells", colour = "Patient")
qc_figure <- (p1 | p2) / (p3 | p4) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(qc_figure, file.path(results_root, "figures", "supplementary", "Figure_S1_QC"), 183, 135)

p_umap <- ggplot(umap, aes(x = umapunintegrated_1, y = umapunintegrated_2, colour = sample_id)) +
  geom_point(size = 0.12, alpha = 0.65) +
  coord_equal() + guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Sample", title = "Unintegrated expression space")
save_pub_r(p_umap, file.path(results_root, "figures", "supplementary", "Figure_S2_unintegrated_UMAP"), 89, 75)

capture_versions(
  c("R", "Matrix", "DropletUtils", "SingleCellExperiment", "scater", "scDblFinder", "SoupX",
    "Seurat", "ggplot2", "patchwork", "svglite", "ragg"),
  file.path(logs_root, "scrna_qc_package_versions.tsv")
)

status <- list(
  stage = "scrna_qc",
  completed = TRUE,
  raw_matrices = nrow(sample_manifest),
  final_cells = sum(qc_summary$final_singlets),
  filter_policy = "EmptyDrops FDR<=0.01; per-sample robust QC on detected genes, UMI counts and mitochondrial fraction; scDblFinder per library",
  ambient_policy = "SoupX per library; failures retained explicitly in ambient_status",
  integration_policy = "No batch integration for the canonical QC object; unintegrated PCA/UMAP only"
)
writeLines(jsonlite::toJSON(status, auto_unbox = TRUE, pretty = TRUE),
           file.path(logs_root, "scrna_qc_status.json"))
message("QC complete: ", sum(qc_summary$final_singlets), " singlets")
