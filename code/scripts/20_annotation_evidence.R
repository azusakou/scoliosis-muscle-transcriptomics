#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(SingleR)
  library(celldex)
  library(yaml)
  library(ggplot2)
  library(patchwork)
})
source("scripts/lib_common.R")
theme_set(theme_nature_contract())

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else "data/processed/scrna/scrna_qc_unintegrated.rds"
processed_root <- if (length(args) >= 2) args[[2]] else "data/processed"
results_root <- if (length(args) >= 3) args[[3]] else "results"

object <- readRDS(input)
DefaultAssay(object) <- "SoupX"
object <- JoinLayers(object, assay = "SoupX")
if (!"data" %in% Layers(object[["SoupX"]])) object <- NormalizeData(object, verbose = FALSE)
Idents(object) <- "seurat_clusters"

message("Computing cluster differential markers for annotation evidence")
markers <- FindAllMarkers(
  object, assay = "SoupX", only.pos = TRUE, min.pct = 0.10,
  logfc.threshold = 0.25, test.use = "wilcox", verbose = FALSE
)
markers$analysis_role <- "annotation_marker_only_not_condition_inference"
write_tsv(markers, file.path(results_root, "tables", "annotation", "cluster_differential_markers.tsv.gz"))

marker_sets <- yaml::read_yaml("config/canonical_markers.yaml")
average <- AverageExpression(object, assays = "SoupX", group.by = "seurat_clusters", layer = "data", verbose = FALSE)$SoupX
colnames(average) <- sub("^g(?=[0-9])", "", colnames(average), perl = TRUE)
all_markers <- unique(unlist(marker_sets))
available <- intersect(all_markers, rownames(average))
scaled <- t(scale(t(as.matrix(average[available, , drop = FALSE]))))
scaled[!is.finite(scaled)] <- 0

score_rows <- list()
for (label in names(marker_sets)) {
  genes <- intersect(marker_sets[[label]], rownames(average))
  for (cluster in colnames(average)) {
    expression <- average[genes, cluster]
    score_rows[[length(score_rows) + 1]] <- data.frame(
      cluster = cluster,
      candidate_label = label,
      marker_score = mean(scaled[genes, cluster], na.rm = TRUE),
      expressed_marker_count = sum(expression > 0.10),
      marker_count_available = length(genes),
      expressed_markers = paste(genes[expression > 0.10], collapse = ","),
      # AverageExpression(layer="data") exponentiates before averaging.
      mean_normalized_expression = mean(expression)
    )
  }
}
score_table <- do.call(rbind, score_rows)
score_table <- score_table[order(as.numeric(score_table$cluster), -score_table$marker_score), ]
write_tsv(score_table, file.path(results_root, "tables", "annotation", "canonical_marker_scores.tsv"))

message("Running cluster-level SingleR reference mapping")
logcounts <- GetAssayData(object, assay = "SoupX", layer = "data")
sce <- SingleCellExperiment(assays = list(logcounts = logcounts))
reference_cache <- file.path(processed_root, "reference_cache", "gypsum")
dir.create(reference_cache, recursive = TRUE, showWarnings = FALSE)
gypsum::cacheDirectory(reference_cache)
reference <- celldex::HumanPrimaryCellAtlasData()
common <- intersect(rownames(sce), rownames(reference))
if (length(common) < 1000) stop("Fewer than 1,000 shared genes with reference")
prediction <- SingleR(
  test = sce[common, ], ref = reference[common, ], labels = reference$label.fine,
  clusters = object$seurat_clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts"
)
reference_table <- data.frame(
  cluster = rownames(prediction),
  reference_label = prediction$labels,
  pruned_reference_label = prediction$pruned.labels,
  delta_next = prediction$delta.next,
  stringsAsFactors = FALSE
)
write_tsv(reference_table, file.path(results_root, "tables", "annotation", "singleR_HPCA_cluster_labels.tsv"))

composition <- as.data.frame(table(
  cluster = object$seurat_clusters,
  sample_id = object$sample_id,
  patient_id = object$patient_id,
  side = object$side
))
composition <- composition[composition$Freq > 0, ]
totals <- aggregate(Freq ~ cluster, composition, sum)
names(totals)[2] <- "cluster_total"
composition <- merge(composition, totals, by = "cluster")
composition$within_cluster_fraction <- composition$Freq / composition$cluster_total
write_tsv(composition, file.path(results_root, "tables", "annotation", "cluster_sample_composition.tsv"))
dominance <- aggregate(within_cluster_fraction ~ cluster, composition, max)
names(dominance)[2] <- "max_single_sample_fraction"
dominance$dominated_by_one_sample <- dominance$max_single_sample_fraction > 0.80
write_tsv(dominance, file.path(results_root, "tables", "annotation", "cluster_sample_dominance.tsv"))

top_score <- do.call(rbind, lapply(split(score_table, score_table$cluster), function(x) x[order(-x$marker_score), ][1:2, ]))
top_marker_rows <- do.call(rbind, lapply(split(markers, markers$cluster), function(x) head(x[order(-x$avg_log2FC), ], 15)))
write_tsv(top_marker_rows, file.path(results_root, "tables", "annotation", "top15_cluster_markers.tsv"))

review <- merge(top_score, reference_table, by = "cluster", all.x = TRUE)
review <- merge(review, dominance, by = "cluster", all.x = TRUE)
write_tsv(review, file.path(results_root, "tables", "annotation", "annotation_review_matrix.tsv"))

saveRDS(object, file.path(processed_root, "scrna", "scrna_annotation_evidence.rds"), compress = FALSE)
message("Annotation evidence complete; manual cluster map remains required")
