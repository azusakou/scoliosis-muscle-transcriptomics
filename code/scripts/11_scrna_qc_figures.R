#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})
source("scripts/lib_common.R")
theme_set(theme_nature_contract())

args <- commandArgs(trailingOnly = TRUE)
revision_replay <- "--revision-replay" %in% args
args <- args[args != "--revision-replay"]
processed_root <- if (length(args) >= 1) args[[1]] else "data/processed"
results_root <- if (length(args) >= 2) args[[2]] else "results"
logs_root <- if (length(args) >= 3) args[[3]] else "logs"
figure_root <- if (revision_replay) "analysis/figures" else
  file.path(results_root, "figures", "supplementary")
figure_prefix <- if (revision_replay) "figure_" else "Figure_"
status_root <- if (revision_replay) file.path(logs_root, "retained_figure_inputs") else logs_root
dir.create(status_root, recursive = TRUE, showWarnings = FALSE)

qc_summary <- utils::read.delim(file.path(results_root, "tables", "qc", "scrna_qc_summary.tsv"))
threshold_table <- utils::read.delim(file.path(results_root, "tables", "qc", "scrna_qc_thresholds.tsv"))
plot_meta <- utils::read.delim(gzfile(file.path(results_root, "tables", "qc", "scrna_cell_qc_metadata.tsv.gz")))
umap <- utils::read.delim(file.path(results_root, "tables", "qc", "scrna_unintegrated_umap.tsv"))
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
save_pub_r(qc_figure, file.path(figure_root, paste0(figure_prefix, "S1_QC")), 183, 135)

p_umap <- ggplot(umap, aes(x = umapunintegrated_1, y = umapunintegrated_2, colour = sample_id)) +
  geom_point(size = 0.12, alpha = 0.65) + coord_equal() +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1))) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Sample", title = "Unintegrated expression space")
save_pub_r(p_umap, file.path(figure_root, paste0(figure_prefix, "S2_unintegrated_UMAP")), 89, 75)

status <- list(
  stage = "scrna_qc_figure_export", completed = TRUE, final_cells = sum(qc_summary$final_singlets),
  scope = "Plots rendered from retained per-cell QC metrics and coordinates; no new droplet calling, SoupX, doublet detection or UMAP fit",
  revision_replay = revision_replay,
  filter_policy = "EmptyDrops FDR<=0.01; per-sample robust QC on detected genes, UMI counts and mitochondrial fraction; scDblFinder per library",
  ambient_policy = "SoupX per library; failures retained explicitly in ambient_status",
  integration_policy = "No batch integration for the canonical QC object; unintegrated PCA/UMAP only",
  figure_backend = "R only"
)
writeLines(jsonlite::toJSON(status, auto_unbox = TRUE, pretty = TRUE), file.path(status_root, "scrna_qc_figure_export.json"))
message("QC figure export complete")
