#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(edgeR); library(limma) })
source("scripts/lib_common.R")
set.seed(324067)

raw_root <- raw_data_root()
raw_dir <- file.path(raw_root, "bulk_gse324067")
meta <- read.delim("metadata/audit/gse324067_samples.tsv", check.names = FALSE,
                   stringsAsFactors = FALSE)
stopifnot(nrow(meta) == 18L,
          all(meta$processed_label_agreement %in% TRUE),
          all(meta$analysis_role == "adjacent_tissue_specificity_negative_control"))

read_count_file <- function(filename) {
  tab <- read.delim(gzfile(file.path(raw_dir, filename)), check.names = FALSE,
                    stringsAsFactors = FALSE)
  stopifnot(ncol(tab) == 2L, names(tab)[1] == "Geneid")
  names(tab) <- c("gene", "count")
  tab
}
tables <- lapply(meta$processed_file, read_count_file)
reference_genes <- tables[[1]]$gene
stopifnot(all(vapply(tables, function(z) identical(z$gene, reference_genes), logical(1))))
counts <- do.call(cbind, lapply(tables, `[[`, "count"))
rownames(counts) <- make.unique(reference_genes)
colnames(counts) <- paste(meta$patient_id, meta$processed_side, sep = "__")
storage.mode(counts) <- "integer"
stopifnot(!anyNA(counts), all(counts >= 0))

meta$sample_id <- colnames(counts)
meta$side <- factor(meta$processed_side, levels = c("non_curve", "convex", "concave"))

# Normalize all three sites together for projection and descriptive comparisons.
keep_all <- filterByExpr(counts, group = meta$side, min.count = 10)
y_all <- calcNormFactors(DGEList(counts[keep_all, , drop = FALSE]))
logcpm_all <- cpm(y_all, log = TRUE, prior.count = 2)
projection_matrix <- label_blind_logcpm(counts, y_all$samples$norm.factors)
projection_logcpm <- projection_matrix$logcpm
pc <- prcomp(t(logcpm_all), scale. = FALSE)
qc <- data.frame(
  sample_id = meta$sample_id, patient_id = meta$patient_id,
  side = as.character(meta$side), processed_file = meta$processed_file,
  geo_accession_attached = meta$geo_accession_attached,
  geo_attachment_conflict = meta$geo_attachment_conflict,
  library_size = colSums(counts), detected_genes = colSums(counts > 0),
  normalization_factor = y_all$samples$norm.factors,
  PC1 = pc$x[meta$sample_id, 1], PC2 = pc$x[meta$sample_id, 2]
)
write_tsv(qc, "results/tables/adjacent_tissue_gse324067/GSE324067_QC.tsv")

# Patient-blocked primary adjacent-tissue contrast: concave minus convex.
paired <- meta$processed_side %in% c("concave", "convex")
meta_pair <- droplevels(meta[paired, ])
meta_pair$side <- factor(meta_pair$processed_side, levels = c("convex", "concave"))
counts_pair <- counts[, paired, drop = FALSE]
keep <- filterByExpr(counts_pair, group = meta_pair$side, min.count = 10)
y <- calcNormFactors(DGEList(counts_pair[keep, , drop = FALSE]))
design <- model.matrix(~factor(patient_id) + side, meta_pair)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "sideconcave")
de <- topTags(qlf, n = Inf, sort.by = "none")$table
de$gene <- rownames(de)
de$contrast <- "osteoblast_concave_minus_convex"
de$n_paired_patients <- length(unique(meta_pair$patient_id))
de$inference_scope <- "adjacent_tissue_specificity_not_muscle_validation"
write_tsv(de, "results/tables/adjacent_tissue_gse324067/GSE324067_osteoblast_paired_edgeR.tsv.gz")

# Within-patient three-site contrasts, retained as contextual specificity analyses.
design3 <- model.matrix(~factor(patient_id) + side, meta)
y3 <- estimateDisp(y_all, design3, robust = TRUE)
fit3 <- glmQLFit(y3, design3, robust = TRUE)
contrast_matrix <- matrix(0, nrow = ncol(design3), ncol = 3,
                          dimnames = list(colnames(design3),
                                          c("concave_minus_noncurve", "convex_minus_noncurve",
                                            "concave_minus_convex")))
contrast_matrix["sideconcave", "concave_minus_noncurve"] <- 1
contrast_matrix["sideconvex", "convex_minus_noncurve"] <- 1
contrast_matrix["sideconcave", "concave_minus_convex"] <- 1
contrast_matrix["sideconvex", "concave_minus_convex"] <- -1
three_site <- do.call(rbind, lapply(colnames(contrast_matrix), function(label) {
  test <- glmQLFTest(fit3, contrast = contrast_matrix[, label])
  tab <- topTags(test, n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab$contrast <- label
  tab
}))
write_tsv(three_site, "results/tables/adjacent_tissue_gse324067/GSE324067_three_site_edgeR.tsv.gz")

# Project the frozen muscle multicellular axis; absence supports tissue specificity,
# while presence indicates a deformity-wide rather than muscle-restricted program.
signature <- read.delim("data/processed/frozen_multicellular_signature.tsv")
projection <- project_frozen_signature(projection_logcpm, signature)
weights <- projection$used_weights
axis_score <- projection$score
usage <- projection$weights
usage$dataset <- "GSE324067_osteoblast"
usage$detectability_rule <- sprintf("CPM>=%.1f_in_at_least_%d_libraries_label_blind",
                                    projection_matrix$min_cpm,
                                    projection_matrix$min_samples)
write_tsv(usage, "results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_projection_genes.tsv")
axis <- data.frame(
  sample_id = names(axis_score), score = as.numeric(axis_score),
  patient_id = meta$patient_id[match(names(axis_score), meta$sample_id)],
  side = meta$processed_side[match(names(axis_score), meta$sample_id)]
)
axis_eff <- merge(axis[axis$side == "concave", ], axis[axis$side == "convex", ],
                  by = "patient_id", suffixes = c(".concave", ".convex"))
axis_eff$paired_difference <- axis_eff$score.concave - axis_eff$score.convex
tt <- t.test(axis_eff$paired_difference)
sign_grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(axis_eff))))
permuted <- as.numeric(sign_grid %*% axis_eff$paired_difference / nrow(axis_eff))
exact_p <- mean(abs(permuted) >= abs(mean(axis_eff$paired_difference)) - 1e-12)
axis_test <- data.frame(
  contrast = "osteoblast_concave_minus_convex",
  n_pairs = nrow(axis_eff), mean_paired_difference = mean(axis_eff$paired_difference),
  se = sd(axis_eff$paired_difference) / sqrt(nrow(axis_eff)),
  t_test_p = tt$p.value, exact_sign_flip_p = exact_p,
  positive_pairs = sum(axis_eff$paired_difference > 0),
  interpretation_scope = "adjacent_tissue_specificity_test"
)
write_tsv(axis, "results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_scores.tsv")
write_tsv(axis_eff, "results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_paired_effects.tsv")
write_tsv(axis_test, "results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_test.tsv")

# Frozen pathway projection with paired effects and exact sign-flip p values.
frozen <- readRDS("data/processed/frozen_pathway_definitions.rds")
path_scores <- gene_set_zscore(logcpm_all, frozen$pathways, min_genes = 5)
path_long <- do.call(rbind, lapply(rownames(path_scores), function(gs) data.frame(
  gene_set = gs, sample_id = colnames(path_scores), score = path_scores[gs, ],
  patient_id = meta$patient_id[match(colnames(path_scores), meta$sample_id)],
  side = meta$processed_side[match(colnames(path_scores), meta$sample_id)]
)))
path_eff <- merge(path_long[path_long$side == "concave", ],
                  path_long[path_long$side == "convex", ],
                  by = c("gene_set", "patient_id"), suffixes = c(".concave", ".convex"))
path_eff$paired_difference <- path_eff$score.concave - path_eff$score.convex
path_tests <- do.call(rbind, lapply(split(path_eff, path_eff$gene_set), function(z) {
  obs <- mean(z$paired_difference)
  grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(z))))
  perm <- as.numeric(grid %*% z$paired_difference / nrow(z))
  data.frame(
    gene_set = z$gene_set[1], n_pairs = nrow(z), mean_paired_difference = obs,
    positive_pairs = sum(z$paired_difference > 0),
    t_test_p = t.test(z$paired_difference)$p.value,
    exact_sign_flip_p = mean(abs(perm) >= abs(obs) - 1e-12)
  )
}))
path_tests$FDR_t_test <- p.adjust(path_tests$t_test_p, "BH")
path_tests$FDR_exact <- p.adjust(path_tests$exact_sign_flip_p, "BH")
write_tsv(path_long, "results/tables/adjacent_tissue_gse324067/frozen_pathway_scores.tsv")
write_tsv(path_eff, "results/tables/adjacent_tissue_gse324067/frozen_pathway_paired_effects.tsv")
write_tsv(path_tests, "results/tables/adjacent_tissue_gse324067/frozen_pathway_paired_tests.tsv")

# Leave-one-patient-out sensitivity for the frozen axis.
loo <- do.call(rbind, lapply(unique(axis_eff$patient_id), function(left_out) {
  z <- axis_eff[axis_eff$patient_id != left_out, ]
  data.frame(left_out_patient = left_out, n_pairs = nrow(z),
             mean_paired_difference = mean(z$paired_difference),
             positive_pairs = sum(z$paired_difference > 0),
             t_test_p = t.test(z$paired_difference)$p.value)
}))
write_tsv(loo, "results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_leave_one_patient_out.tsv")

# Cross-tissue effect concordance is descriptive and is never meta-analysed with muscle.
muscle <- read.delim(gzfile(analysis_input("results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz")))
cross <- merge(muscle[, c("gene", "logFC", "PValue", "FDR")],
               de[, c("gene", "logFC", "PValue", "FDR")], by = "gene",
               suffixes = c(".muscle_GSE254300", ".osteoblast_GSE324067"))
cross$sign_concordant <- sign(cross$logFC.muscle_GSE254300) == sign(cross$logFC.osteoblast_GSE324067)
summary <- data.frame(
  metric = c("spearman_logFC_all", "sign_concordance_all", "sign_concordance_muscle_FDR05"),
  value = c(
    cor(cross$logFC.muscle_GSE254300, cross$logFC.osteoblast_GSE324067, method = "spearman"),
    mean(cross$sign_concordant),
    mean(cross$sign_concordant[cross$FDR.muscle_GSE254300 < 0.05])
  ),
  n = c(nrow(cross), nrow(cross), sum(cross$FDR.muscle_GSE254300 < 0.05)),
  interpretation = "descriptive_cross_tissue_specificity_not_replication"
)
write_tsv(cross, "results/tables/adjacent_tissue_gse324067/muscle_osteoblast_effect_comparison.tsv.gz")
write_tsv(summary, "results/tables/adjacent_tissue_gse324067/muscle_osteoblast_effect_summary.tsv")

capture_versions(c("edgeR", "limma"), "logs/gse324067_package_versions.tsv")
message("GSE324067 adjacent-tissue analysis complete for six patient triplets.")
