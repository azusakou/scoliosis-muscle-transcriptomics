#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("scripts/lib_common.R")

# Figure contract
# Core conclusion: original source-labelled muscle contrasts are exploratory;
# adjacent tissue and biopsy-to-culture comparisons have separate roles.
# Archetype: quantitative grid; patient-level data and explicit RNA-gate flags.
# Output: double-column 183 x 130 mm; editable SVG/PDF and 600-dpi TIFF.

revision_replay <- "--revision-replay" %in% commandArgs(trailingOnly=TRUE)
out <- if (revision_replay) "analysis/figures/figure_S3_dataset_extension" else
  "results/figures/supplementary/Figure_S3_dataset_extension"
source_dir <- "analysis/tables"
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

axis12 <- read.delim("results/tables/meta_validation_requant/frozen_axis_12_patient_effects.tsv")
axis12_test <- read.delim("results/tables/meta_validation_requant/frozen_axis_12_patient_test.tsv")
axis12$cohort_label <- ifelse(axis12$cohort == "GSE254300", "GSE254300 (n=5)",
                              "GSE311180 raw re-quant. (n=7)")
rna_pairs <- read.delim(analysis_input("analysis/audit_gse311180_rna_expected_pairs.tsv"))
failed <- rna_pairs$patient_id[rna_pairs$interpretation == "does_not_meet_consistency_gates"]
stopifnot(setequal(failed, c("STB-084", "STB-100")))
axis12$rna_gate_status <- ifelse(axis12$cohort == "GSE254300", "not_assessed_here",
                               ifelse(axis12$patient_id %in% failed, "nonpassing", "passing"))
axis12$patient_plot <- paste0(axis12$patient_id, ifelse(axis12$rna_gate_status == "nonpassing", "*", ""))
axis12$direction <- axis12$paired_difference > 0
p_a <- ggplot(axis12, aes(reorder(patient_plot, paired_difference), paired_difference,
                          fill = direction)) +
  geom_col(width = .7) + geom_hline(yintercept = 0, linewidth = .35) +
  facet_grid(~cohort_label, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(`TRUE` = palette_contract[["concave"]],
                               `FALSE` = palette_contract[["convex"]]), guide = "none") +
  labs(x = "Patient", y = "Signed-score paired difference",
       title = "Original source-label score analysis",
       subtitle = sprintf("5 + 7 pairs; %d/%d positive; exact P = %.4f",
                          axis12_test$positive_pairs, axis12_test$n_patients,
                          axis12_test$exact_sign_flip_p),
       caption = "* Two RNA-gate failures included; exclusion removes score significance.\nUnits: mean signed gene-wise z-score, concave minus convex.") +
  theme_nature_contract(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(), strip.placement = "outside")

replication <- read.delim(gzfile("results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz"))
rep_summary <- read.delim("results/tables/meta_validation_requant/gene_replication_summary.tsv")
rho <- rep_summary$value[rep_summary$metric == "spearman_logFC_all"]
p_b <- ggplot(replication, aes(logFC.GSE254300, logFC.GSE311180_requant)) +
  geom_hline(yintercept = 0, colour = "#D8D8D8", linewidth = .3) +
  geom_vline(xintercept = 0, colour = "#D8D8D8", linewidth = .3) +
  geom_point(size = .55, alpha = .35, colour = "#737D88") +
  labs(x = "GSE254300 fitted log2 fold change", y = "GSE311180 fitted log2 fold change",
       title = "Original all-label gene contrasts",
       subtitle = sprintf("12,995 common genes; Spearman rho = %.3f", rho),
       caption = "All genes shown; no Stouffer-significance selection.\nRNA-gate-passing sensitivity is shown in Figure S10.") +
  theme_nature_contract(7)

osteo <- read.delim("results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_paired_effects.tsv")
osteo_test <- read.delim("results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_test.tsv")
p_c <- ggplot(osteo, aes(reorder(patient_id, paired_difference), paired_difference,
                         fill = paired_difference > 0)) +
  geom_col(width = .68) + geom_hline(yintercept = 0, linewidth = .35) +
  scale_fill_manual(values = c(`TRUE` = "#D9A39C", `FALSE` = "#9ABAD0"), guide = "none") +
  labs(x = "Patient", y = "Signed-score paired difference",
       title = "Adjacent osteoblast-culture comparison",
       subtitle = sprintf("%d/%d positive; exact P = %.3f",
                          osteo_test$positive_pairs, osteo_test$n_pairs,
                          osteo_test$exact_sign_flip_p),
       caption = "Different tissue/culture system; not bilateral muscle validation.") +
  theme_nature_contract(7)

reference <- read.delim("results/tables/reference_gse17503/cell_type_marker_paired_tests.tsv")
reference$label <- gsub(" cell$", "", reference$cell_type)
p_d <- ggplot(reference, aes(mean_biopsy_minus_culture,
                             reorder(label, mean_biopsy_minus_culture))) +
  geom_vline(xintercept = 0, linewidth = .35) +
  geom_segment(aes(x = 0, xend = mean_biopsy_minus_culture, yend = label),
               colour = "#B8B8B8", linewidth = .45) +
  geom_point(aes(fill = positive_pairs == n_pairs), shape = 21, size = 2,
             colour = "#272727", stroke = .3) +
  scale_fill_manual(values = c(`TRUE` = palette_contract[["signal_teal"]],
                               `FALSE` = "white"), guide = "none") +
  labs(x = "Biopsy - matched culture marker score", y = NULL,
       title = "Biopsy-to-culture reference sensitivity",
       subtitle = "Filled: positive biopsy-minus-culture change in all five pairs",
       caption = "Marker-score units; no deposited concave/convex labels.") +
  theme_nature_contract(7)

figure <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(1.1, 1), widths = c(1.08, .92)) +
  plot_annotation(
    title = "Original muscle contrasts and role-restricted dataset context",
    tag_levels = "a"
  ) & theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(figure, out, 183, 158)
ragg::agg_png(paste0(out, ".png"), width=183, height=158, units="mm", res=180)
print(figure); grDevices::dev.off()

write_tsv(axis12, file.path(source_dir, "figure_S3_original_source_label_patient_effects.tsv"))
write_tsv(replication[,c("gene", "logFC.GSE254300", "logFC.GSE311180_requant")],
          file.path(source_dir, "figure_S3_original_common_gene_effects.tsv"))
write_tsv(osteo, file.path(source_dir, "figure_S3_osteoblast_paired_effects.tsv"))
write_tsv(reference[,c("cell_type", "n_genes", "n_pairs", "mean_biopsy_minus_culture", "positive_pairs")],
          file.path(source_dir, "figure_S3_biopsy_culture_marker_effects.tsv"))
message("Figure S3 extension evidence exported using R only.")
