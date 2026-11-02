#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
})
source("scripts/lib_common.R")

set.seed(722100)
object <- readRDS("data/processed/scrna/scrna_annotated.rds")
object <- JoinLayers(object, assay = "SoupX")
counts <- LayerData(object, assay = "SoupX", layer = "counts")
meta <- object@meta.data
meta$pb_id <- paste(meta$patient_id, meta$side, sep = "__")

cell_counts <- as.data.frame(table(meta$cell_type, meta$pb_id), stringsAsFactors = FALSE)
names(cell_counts) <- c("cell_type", "pb_id", "n_cells")
write_tsv(cell_counts, "results/tables/pseudobulk/cell_counts_per_pseudobulk.tsv")

min_cells <- 20L
cell_types <- sort(unique(meta$cell_type))
de_all <- list(); donor_all <- list(); expression_all <- list()

for (ct in cell_types) {
  md <- meta[meta$cell_type == ct, , drop = FALSE]
  tab <- table(md$pb_id)
  valid <- names(tab)[tab >= min_cells]
  if (length(valid) != 4L) next
  pb <- vapply(valid, function(id) Matrix::rowSums(counts[, rownames(md)[md$pb_id == id], drop = FALSE]),
               numeric(nrow(counts)))
  rownames(pb) <- rownames(counts)
  sample_info <- data.frame(
    pb_id = colnames(pb),
    patient_id = sub("__.*$", "", colnames(pb)),
    side = sub("^.*__", "", colnames(pb)), stringsAsFactors = FALSE
  )
  sample_info$side <- factor(sample_info$side, levels = c("convex", "concave"))
  design <- model.matrix(~factor(patient_id) + side, sample_info)
  keep <- filterByExpr(pb, design = design, min.count = 5)
  y <- DGEList(pb[keep, , drop = FALSE])
  y <- calcNormFactors(y)
  stopifnot(qr(design)$rank == ncol(design))
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "sideconcave")
  de <- topTags(qlf, n = Inf, sort.by = "none")$table
  de$gene <- rownames(de)
  de$cell_type <- ct
  de$n_donors <- 2L
  de$contrast <- "concave_minus_convex"
  de$count_assay <- "SoupX_ambient_corrected_counts"
  de$inferential_status <- "exploratory_effect_map_n2_not_population_inference"

  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  donor_effects <- do.call(rbind, lapply(unique(sample_info$patient_id), function(pid) {
    cc <- sample_info$pb_id[sample_info$patient_id == pid & sample_info$side == "concave"]
    cv <- sample_info$pb_id[sample_info$patient_id == pid & sample_info$side == "convex"]
    data.frame(gene = rownames(logcpm), cell_type = ct, patient_id = pid,
               paired_logCPM_difference = logcpm[, cc] - logcpm[, cv])
  }))
  wide <- reshape(donor_effects, idvar = c("gene", "cell_type"), timevar = "patient_id", direction = "wide")
  effect_cols <- grep("^paired_logCPM_difference", names(wide), value = TRUE)
  wide$donor_sign_concordant <- apply(wide[, effect_cols, drop = FALSE], 1, function(z) all(z > 0) || all(z < 0))
  wide$mean_donor_effect <- rowMeans(wide[, effect_cols, drop = FALSE])
  de <- merge(de, wide, by = c("gene", "cell_type"), all.x = TRUE)

  de_all[[ct]] <- de
  donor_all[[ct]] <- donor_effects
  expression_all[[ct]] <- list(logcpm = logcpm, sample_info = sample_info)
  write_tsv(data.frame(gene = rownames(pb), pb),
            file.path("results/tables/pseudobulk/counts", paste0(safe_name(ct), ".tsv.gz")))
}

de_all <- do.call(rbind, de_all)
donor_all <- do.call(rbind, donor_all)
write_tsv(de_all, "results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz")
write_tsv(donor_all, "results/tables/pseudobulk/scrna_donor_specific_effects.tsv.gz")
saveRDS(expression_all, "data/processed/scrna/pseudobulk_logcpm_by_cell_type.rds")

# Patient-level composition effects. With n=2, these are descriptive; no cell-level test is run.
composition <- read.delim("results/tables/composition/cell_type_counts_and_proportions.tsv", check.names = FALSE)
composition_wide <- reshape(composition[, c("patient_id", "side", "cell_type", "proportion")],
                            idvar = c("patient_id", "cell_type"), timevar = "side", direction = "wide")
composition_wide$concave_minus_convex <- composition_wide$proportion.concave - composition_wide$proportion.convex
composition_wide$inferential_status <- "descriptive_n2_no_cell_level_test"
write_tsv(composition_wide, "results/tables/composition/paired_composition_effects.tsv")

capture_versions(c("Seurat", "Matrix", "edgeR", "limma"), "logs/pseudobulk_package_versions.tsv")
message("Patient-level SoupX pseudobulk complete: ", length(expression_all),
        " cell types with >=20 cells in all four samples.")
