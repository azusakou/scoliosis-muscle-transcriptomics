#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
})
source("scripts/lib_common.R")

in_dir <- "results/tables/disease_control_gse254300_cs"
revision_replay <- "--revision-replay" %in% commandArgs(trailingOnly=TRUE)
fig_dir <- if (revision_replay) "analysis/figures" else "results/figures/supplementary"
figure_stem <- if (revision_replay) "figure_S4_GSE254300_CS_specificity" else
  "Figure_S4_GSE254300_CS_specificity"
source_dir <- "analysis/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

# Figure contract
# Core conclusion: the diagnosis comparison is inconclusive; neither diagnosis
# has a detectable score effect, and wide intervals do not establish sharing.
# Archetype: asymmetric quantitative comparison figure.
# Evidence: patient pairing, descriptive gene contrasts, model uncertainty.
# Backend/export: R only; 183-mm SVG/PDF plus 600-dpi TIFF.

COL_AIS <- "#355C8A"
COL_CS <- "#C56A48"
COL_DARK <- "#252A31"
COL_MID <- "#858B94"
COL_LIGHT <- "#D9DDE2"
COL_ACCENT <- "#B23A48"

theme_topbio <- function(base_size = 7.3) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.32, colour = COL_DARK),
      axis.ticks = element_line(linewidth = 0.32, colour = COL_DARK),
      axis.ticks.length = unit(1.5, "mm"),
      axis.title = element_text(size = base_size, colour = COL_DARK),
      axis.text = element_text(size = base_size - 0.5, colour = COL_DARK),
      plot.title = element_text(size = base_size + 0.5, face = "bold",
                                colour = COL_DARK, margin = margin(b = 2, l = 11)),
      plot.subtitle = element_text(size = base_size - 0.15, colour = "#555B63",
                                   lineheight = 1.05, margin = margin(b = 5)),
      plot.caption = element_text(size = base_size - 0.65, colour = "#555B63",
                                  hjust = 0, margin = margin(t = 4)),
      plot.tag = element_text(size = 8.2, face = "bold", colour = COL_DARK),
      plot.tag.position = c(0, 1),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 0.1, face = "bold", colour = COL_DARK),
      legend.title = element_text(size = base_size - 0.5),
      legend.text = element_text(size = base_size - 0.8),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(3.2, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(6, 7, 5, 7)
    )
}

# a, hero panel: preserve every within-patient pair and place cohort means over
# the patient trajectories. Side and mean symbols have an explicit legend.
axis <- read.delim(file.path(in_dir, "frozen_axis_joint_scores.tsv"))
axis_tests <- read.delim(file.path(in_dir, "frozen_axis_within_cohort_tests.tsv"))
axis_int <- read.delim(file.path(in_dir, "frozen_axis_AIS_CS_interaction_test.tsv"))
axis$cohort <- factor(axis$cohort, levels = c("AIS", "CS"))
axis$side <- factor(axis$side, levels = c("convex", "concave"))
axis$cohort_colour <- ifelse(axis$cohort == "AIS", COL_AIS, COL_CS)
axis_mean <- aggregate(score ~ cohort + side, axis, mean)

p_a <- ggplot(axis, aes(side, score, group = pair_id)) +
  geom_hline(yintercept = 0, colour = COL_LIGHT, linewidth = 0.35) +
  geom_line(aes(colour = cohort), linewidth = 0.55, alpha = 0.42) +
  geom_point(aes(colour = cohort, shape = side), size = 2.05, stroke = 0.65) +
  geom_line(data = axis_mean, aes(side, score, group = cohort, colour = cohort),
            linewidth = 1.05, inherit.aes = FALSE) +
  geom_point(data = axis_mean, aes(side, score, colour = cohort, shape = "Mean"),
             size = 3.0, inherit.aes = FALSE) +
  facet_grid(cohort ~ ., switch = "y") +
  scale_colour_manual(values = c(AIS = COL_AIS, CS = COL_CS), guide = "none") +
  scale_fill_manual(values = c(AIS = COL_AIS, CS = COL_CS), guide = "none") +
  scale_shape_manual(values = c(convex = 1, concave = 16, Mean = 18), name = NULL,
                     labels = c(convex = "Convex", concave = "Concave", Mean = "Mean")) +
  scale_x_discrete(labels = c(convex = "Convex", concave = "Concave")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.11))) +
  labs(
    tag = "a",
    x = NULL, y = "Mean signed gene-wise z-score",
    title = "Diagnosis comparison remains unresolved",
    subtitle = sprintf("AIS: %d/5 positive, P = %.4f\nCS: %d/5 positive, P = %.4f",
                       axis_tests$positive_pairs[axis_tests$cohort == "AIS"],
                       axis_tests$exact_sign_flip_p[axis_tests$cohort == "AIS"],
                       axis_tests$positive_pairs[axis_tests$cohort == "CS"],
                       axis_tests$exact_sign_flip_p[axis_tests$cohort == "CS"])
  ) +
  theme_topbio(7.5) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
        panel.spacing.y = unit(4.5, "mm"), legend.position = "bottom"
  )

# b, genome-wide density. The binned background avoids an opaque point cloud;
# named genes are overlaid as exact effects, and y=x shows identical effects.
gene <- read.delim(gzfile(file.path(in_dir, "GSE254300_AIS_CS_gene_effect_comparison.tsv.gz")))
gene_summary <- read.delim(file.path(in_dir, "GSE254300_AIS_CS_gene_effect_summary.tsv"))
rho <- gene_summary$value[gene_summary$metric == "spearman_AIS_CS_logFC"]
concordance <- gene_summary$value[gene_summary$metric == "sign_concordance_all"]
gene$highlight <- gene$gene %in% c("ANGPT1", "LEP")
lim <- max(abs(c(gene$logFC.AIS, gene$logFC.CS)), na.rm = TRUE)
lim <- ceiling(lim * 2) / 2

p_b <- ggplot(gene, aes(logFC.AIS, logFC.CS)) +
  geom_hline(yintercept = 0, colour = COL_LIGHT, linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = COL_LIGHT, linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = COL_MID, linewidth = 0.4) +
  geom_bin_2d(bins = 48, colour = NA) +
  scale_fill_gradient(low = "#EEF1F4", high = "#294B73", trans = "sqrt",
                      guide = "none") +
  geom_point(data = gene[gene$highlight, ], colour = COL_ACCENT,
             fill = "white", shape = 21, size = 2.4, stroke = 0.8) +
  geom_text_repel(data = gene[gene$highlight, ], aes(label = gene),
                  colour = COL_DARK, size = 2.45, fontface = "bold",
                  min.segment.length = 0, segment.colour = COL_MID,
                  segment.size = 0.3, box.padding = 0.35, point.padding = 0.25,
                  seed = 254300, max.overlaps = Inf) +
  annotate("text", x = -lim * 0.94, y = lim * 0.92,
           label = sprintf("rho = %.3f\nsign concordance = %.1f%%",
                           rho, 100 * concordance),
           hjust = 0, vjust = 1, size = 2.25, colour = COL_DARK, lineheight = 1.10) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
  labs(tag = "b", x = "AIS paired log2 fold-change", y = "CS paired log2 fold-change",
       title = "Descriptive diagnosis-level gene contrasts",
       subtitle = "Dashed line: y = x; no equivalence is inferred") +
  theme_topbio(7.2)

# c, one aligned forest panel: model sensitivity within each diagnosis. The
# diagnosis interaction is reported as an exact patient-label test in caption.
models <- read.delim(file.path(in_dir, "frozen_axis_model_sensitivity.tsv"))
stopifnot(nrow(models) == 6L, all(is.finite(models$residual_df)), all(models$residual_df > 0))
models$model_id <- models$model
models$lower <- models$estimate_concave_minus_convex - qt(.975, models$residual_df) * models$se
models$upper <- models$estimate_concave_minus_convex + qt(.975, models$residual_df) * models$se
models$model <- factor(
  models$model,
  levels = c("biology_plus_composition_PC1", "biology_plus_library_QC", "biology_only"),
  labels = c("+ composition PC1", "+ library size", "Patient + side")
)
models$cohort <- factor(models$cohort, levels = c("AIS", "CS"))

p_c <- ggplot(models, aes(estimate_concave_minus_convex, model, colour = cohort)) +
  geom_vline(xintercept = 0, colour = COL_DARK, linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                orientation = "y", width = 0, linewidth = 0.55,
                position = position_dodge(width = 0.34)) +
  geom_point(size = 2.25, position = position_dodge(width = 0.34)) +
  scale_colour_manual(values = c(AIS = COL_AIS, CS = COL_CS), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(
    tag = "c",
    x = "Adjusted paired effect (95% t confidence interval)", y = NULL,
    title = "Wide intervals limit interpretation",
    subtitle = "Student-t quantiles use model residual df (3 or 4)",
    caption = sprintf("AIS-CS difference = %.3f; exact P = %.3f.\nFailure to detect an interaction\ndoes not establish equivalence.",
                      axis_int$difference_in_means,
                      axis_int$exact_patient_label_permutation_p)
  ) +
  theme_topbio(7.2) +
  theme(legend.position = "top", legend.justification = "left",
        legend.box.margin = margin(0, 0, -2, 0))

# Asymmetric page: the patient-level panel spans both rows; the two narrower
# panels display descriptive gene contrasts and model uncertainty.
figure <- p_a | (p_b / p_c + plot_layout(heights = c(1.32, 0.68))) +
  plot_layout(widths = c(1.05, 1.25)) +
  plot_annotation() &
  theme(plot.tag = element_text(size = 8.2, face = "bold", colour = COL_DARK),
        plot.tag.position = c(0.01, 0.99))

save_pub_r(figure, file.path(fig_dir, figure_stem), 183, 145)
ragg::agg_png(paste0(file.path(fig_dir, figure_stem), ".png"),
              width=183, height=145, units="mm", res=180)
print(figure); grDevices::dev.off()

write_tsv(axis, file.path(source_dir, "figure_S4_joint_signed_scores.tsv"))
write_tsv(axis_mean, file.path(source_dir, "figure_S4_joint_side_means.tsv"))
write_tsv(gene[,c("gene", "logFC.AIS", "logFC.CS")],
          file.path(source_dir, "figure_S4_descriptive_gene_effects.tsv"))
write_tsv(models[,c("cohort", "model_id", "estimate_concave_minus_convex", "se",
                   "residual_df", "lower", "upper")],
          file.path(source_dir, "figure_S4_model_t_confidence_intervals.tsv"))
message("Supplementary Figure S4 exported in top-journal asymmetric layout using R only.")
