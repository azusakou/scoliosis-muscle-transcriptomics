#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readxl)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(Seurat)
  library(nnls)
})
source("scripts/lib_common.R")
set.seed(254305)

out_dir <- "results/tables/disease_control_gse254300_cs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_root <- raw_data_root()
bulk_file <- file.path(raw_root, "bulk", "GSE254300_readcount.xlsx")
audit <- read.delim("metadata/audit/bulk_samples.tsv", check.names = FALSE)
audit <- audit[audit$cohort %in% c("AIS", "congenital_scoliosis"), ]
stopifnot(nrow(audit) == 20L,
          all(table(audit$cohort) == 10L),
          all(table(audit$cohort, audit$side) == 5L),
          all(table(audit$patient_id) == 2L))

raw <- as.data.frame(read_excel(bulk_file))
stopifnot(all(audit$processed_sample_id %in% names(raw)))
raw$GeneID <- make.unique(as.character(raw$GeneID))
counts_all <- as.matrix(raw[, audit$processed_sample_id, drop = FALSE])
rownames(counts_all) <- raw$GeneID
storage.mode(counts_all) <- "numeric"
info_all <- audit[match(colnames(counts_all), audit$processed_sample_id), ]
info_all$cohort_label <- ifelse(info_all$cohort == "AIS", "AIS", "CS")
info_all$pair_id <- paste(info_all$cohort_label, info_all$patient_id, sep = "_")
info_all$side <- factor(info_all$side, levels = c("convex", "concave"))

# Joint TMM normalization is used only for the AIS-versus-CS specificity analysis.
# The original AIS validation outputs remain unchanged.
keep_all <- filterByExpr(counts_all, group = interaction(info_all$cohort_label, info_all$side), min.count = 10)
y_all <- DGEList(counts_all[keep_all, , drop = FALSE])
y_all <- calcNormFactors(y_all)
logcpm_all <- cpm(y_all, log = TRUE, prior.count = 2)
projection_matrix <- label_blind_logcpm(counts_all, y_all$samples$norm.factors)
projection_logcpm <- projection_matrix$logcpm

qc <- data.frame(
  sample_id = colnames(counts_all), patient_id = info_all$patient_id,
  pair_id = info_all$pair_id, cohort = info_all$cohort_label,
  side = as.character(info_all$side), library_size = colSums(counts_all),
  detected_genes = colSums(counts_all > 0),
  normalization_factor = y_all$samples$norm.factors
)
pc <- prcomp(t(logcpm_all), scale. = FALSE)
qc$PC1 <- pc$x[qc$sample_id, 1]
qc$PC2 <- pc$x[qc$sample_id, 2]
write_tsv(qc, file.path(out_dir, "GSE254300_AIS_CS_joint_QC.tsv"))
write_tsv(data.frame(gene = rownames(logcpm_all), logcpm_all),
          file.path(out_dir, "GSE254300_AIS_CS_joint_TMM_logCPM.tsv.gz"))

# Paired CS differential expression.
cs_idx <- info_all$cohort_label == "CS"
counts_cs <- counts_all[, cs_idx, drop = FALSE]
info_cs <- droplevels(info_all[cs_idx, ])
keep_cs <- filterByExpr(counts_cs, group = info_cs$side, min.count = 10)
y_cs <- DGEList(counts_cs[keep_cs, , drop = FALSE])
y_cs <- calcNormFactors(y_cs)
design_cs <- model.matrix(~factor(patient_id) + side, info_cs)
stopifnot(qr(design_cs)$rank == ncol(design_cs))
y_cs <- estimateDisp(y_cs, design_cs, robust = TRUE)
fit_cs <- glmQLFit(y_cs, design_cs, robust = TRUE)
qlf_cs <- glmQLFTest(fit_cs, coef = "sideconcave")
de_cs <- topTags(qlf_cs, n = Inf, sort.by = "none")$table
de_cs$gene <- rownames(de_cs)
de_cs$contrast <- "CS_concave_minus_convex"
de_cs$n_paired_patients <- 5L
write_tsv(de_cs, file.path(out_dir, "GSE254300_CS_paired_edgeR_full.tsv.gz"))

# Gene-level AIS-by-side interaction. Patient fixed effects absorb cohort main effects;
# the explicit interaction tests whether paired side effects differ between AIS and CS.
info_all$side_concave <- as.integer(info_all$side == "concave")
info_all$cohort_AIS <- as.integer(info_all$cohort_label == "AIS")
design_int <- model.matrix(~factor(pair_id) + side_concave + side_concave:cohort_AIS, info_all)
stopifnot(qr(design_int)$rank == ncol(design_int))
y_int <- estimateDisp(y_all, design_int, robust = TRUE)
fit_int <- glmQLFit(y_int, design_int, robust = TRUE)
interaction_coef <- grep("side_concave:cohort_AIS", colnames(design_int), fixed = TRUE)
stopifnot(length(interaction_coef) == 1L)
qlf_int <- glmQLFTest(fit_int, coef = interaction_coef)
de_int <- topTags(qlf_int, n = Inf, sort.by = "none")$table
de_int$gene <- rownames(de_int)
de_int$contrast <- "AIS_minus_CS_difference_in_concave_minus_convex_effect"
de_int$n_AIS_pairs <- 5L
de_int$n_CS_pairs <- 5L
write_tsv(de_int, file.path(out_dir, "GSE254300_AIS_by_side_interaction_edgeR_full.tsv.gz"))

# Cross-diagnosis gene-effect comparison on the common tested gene set.
de_ais <- read.delim(gzfile(analysis_input("results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz")))
gene_compare <- merge(
  de_ais[, c("gene", "logFC", "PValue", "FDR")],
  de_cs[, c("gene", "logFC", "PValue", "FDR")], by = "gene",
  suffixes = c(".AIS", ".CS")
)
gene_compare <- merge(gene_compare, de_int[, c("gene", "logFC", "PValue", "FDR")],
                      by = "gene", all.x = TRUE)
names(gene_compare)[names(gene_compare) == "logFC"] <- "interaction_logFC_AIS_minus_CS"
names(gene_compare)[names(gene_compare) == "PValue"] <- "interaction_PValue"
names(gene_compare)[names(gene_compare) == "FDR"] <- "interaction_FDR"
gene_compare$sign_concordant <- sign(gene_compare$logFC.AIS) == sign(gene_compare$logFC.CS)
gene_summary <- data.frame(
  metric = c("common_tested_genes", "spearman_AIS_CS_logFC", "sign_concordance_all",
             "CS_genes_FDR05_common", "interaction_genes_FDR05"),
  value = c(nrow(gene_compare),
            cor(gene_compare$logFC.AIS, gene_compare$logFC.CS, method = "spearman"),
            mean(gene_compare$sign_concordant, na.rm = TRUE),
            sum(gene_compare$FDR.CS < 0.05, na.rm = TRUE),
            sum(gene_compare$interaction_FDR < 0.05, na.rm = TRUE))
)
write_tsv(gene_compare, file.path(out_dir, "GSE254300_AIS_CS_gene_effect_comparison.tsv.gz"))
write_tsv(gene_summary, file.path(out_dir, "GSE254300_AIS_CS_gene_effect_summary.tsv"))

# Frozen discovery-only signature projection in the jointly normalized matrix.
signature <- read.delim("data/processed/frozen_multicellular_signature.tsv")
projection <- project_frozen_signature(projection_logcpm, signature)
gene_weights <- projection$used_weights
axis_score <- projection$score
usage <- projection$weights
usage$dataset <- "GSE254300_AIS_CS_joint"
usage$detectability_rule <- sprintf("CPM>=%.1f_in_at_least_%d_libraries_label_blind",
                                    projection_matrix$min_cpm,
                                    projection_matrix$min_samples)
write_tsv(usage, file.path(out_dir, "frozen_axis_projection_genes.tsv"))
axis <- data.frame(
  sample_id = names(axis_score), score = as.numeric(axis_score),
  patient_id = info_all$patient_id[match(names(axis_score), info_all$processed_sample_id)],
  pair_id = info_all$pair_id[match(names(axis_score), info_all$processed_sample_id)],
  cohort = info_all$cohort_label[match(names(axis_score), info_all$processed_sample_id)],
  side = as.character(info_all$side[match(names(axis_score), info_all$processed_sample_id)]),
  n_signature_genes = nrow(gene_weights)
)
axis_eff <- merge(axis[axis$side == "concave", ], axis[axis$side == "convex", ],
                  by = c("pair_id", "patient_id", "cohort"), suffixes = c(".concave", ".convex"))
axis_eff$paired_difference <- axis_eff$score.concave - axis_eff$score.convex

exact_sign_flip <- function(x) {
  n <- length(x)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  null <- as.vector(signs %*% x / n)
  observed <- mean(x)
  data.frame(n_pairs = n, positive_pairs = sum(x > 0), mean_paired_difference = observed,
             median_paired_difference = median(x), paired_t_p = t.test(x)$p.value,
             exact_sign_flip_p = mean(abs(null) >= abs(observed) - 1e-12),
             min_attainable_two_sided_p = 2 / (2^n))
}
axis_tests <- do.call(rbind, lapply(split(axis_eff, axis_eff$cohort), function(z) {
  cbind(cohort = z$cohort[1], exact_sign_flip(z$paired_difference))
}))

# Exact diagnosis-label permutation on patient-level paired differences.
observed_interaction <- with(axis_eff,
  mean(paired_difference[cohort == "AIS"]) - mean(paired_difference[cohort == "CS"]))
ais_combinations <- combn(seq_len(nrow(axis_eff)), 5L)
null_interaction <- apply(ais_combinations, 2, function(idx) {
  mean(axis_eff$paired_difference[idx]) - mean(axis_eff$paired_difference[-idx])
})
axis_interaction <- data.frame(
  contrast = "AIS_minus_CS_paired_axis_difference", n_AIS_pairs = 5L, n_CS_pairs = 5L,
  mean_AIS = mean(axis_eff$paired_difference[axis_eff$cohort == "AIS"]),
  mean_CS = mean(axis_eff$paired_difference[axis_eff$cohort == "CS"]),
  difference_in_means = observed_interaction,
  exact_patient_label_permutation_p = mean(abs(null_interaction) >= abs(observed_interaction) - 1e-12),
  n_label_permutations = length(null_interaction)
)
write_tsv(axis, file.path(out_dir, "frozen_axis_joint_scores.tsv"))
write_tsv(axis_eff, file.path(out_dir, "frozen_axis_joint_paired_effects.tsv"))
write_tsv(axis_tests, file.path(out_dir, "frozen_axis_within_cohort_tests.tsv"))
write_tsv(axis_interaction, file.path(out_dir, "frozen_axis_AIS_CS_interaction_test.tsv"))
write_tsv(data.frame(permuted_difference = null_interaction),
          file.path(out_dir, "frozen_axis_AIS_CS_interaction_null.tsv"))

# Composition sensitivity using the same frozen scRNA marker-reference recipe
# as the primary AIS bulk analysis. Estimates are sensitivity covariates, not
# absolute histological fractions.
obj <- readRDS("data/processed/scrna/scrna_annotated.rds")
obj <- JoinLayers(obj, assay = "SoupX")
sc <- LayerData(obj, assay = "SoupX", layer = "counts")
ct <- as.character(obj$cell_type)
eligible <- names(which(table(ct) >= 50))
reference <- sapply(eligible, function(k) Matrix::rowSums(sc[, ct == k, drop = FALSE]))
reference <- sweep(reference, 2, colSums(reference), "/") * 1e6
reference_log <- log2(reference + 1)
markers <- read.delim("results/tables/bulk/deconvolution_marker_reference.tsv")
marker_genes <- intersect(unique(markers$gene), intersect(rownames(reference_log), rownames(logcpm_all)))
ref_use <- reference_log[marker_genes, eligible, drop = FALSE]
bulk_use <- logcpm_all[marker_genes, , drop = FALSE]
deconv_one <- function(b, method) {
  if (method == "NNLS") {
    x <- nnls::nnls(ref_use, b)$x
  } else {
    lambda <- 1
    fn <- function(x) sum((b - ref_use %*% x)^2) + lambda * sum(x^2)
    x <- optim(rep(1 / ncol(ref_use), ncol(ref_use)), fn,
               method = "L-BFGS-B", lower = 0)$par
  }
  if (sum(x) > 0) x / sum(x) else x
}
props <- do.call(rbind, lapply(c("NNLS", "ridge_NNLS"), function(method) {
  do.call(rbind, lapply(colnames(bulk_use), function(s) {
    x <- deconv_one(bulk_use[, s], method)
    data.frame(method = method, sample_id = s, cell_type = eligible, proportion = x)
  }))
}))
props$patient_id <- info_all$patient_id[match(props$sample_id, info_all$processed_sample_id)]
props$pair_id <- info_all$pair_id[match(props$sample_id, info_all$processed_sample_id)]
props$cohort <- info_all$cohort_label[match(props$sample_id, info_all$processed_sample_id)]
props$side <- as.character(info_all$side[match(props$sample_id, info_all$processed_sample_id)])
write_tsv(props, file.path(out_dir, "deconvolution_two_methods.tsv"))
wide_props <- reshape(props,
                      idvar = c("sample_id", "pair_id", "patient_id", "cohort", "side", "cell_type"),
                      timevar = "method", direction = "wide")
agreement <- do.call(rbind, lapply(split(wide_props, wide_props$cell_type), function(z) {
  data.frame(cell_type = z$cell_type[1],
             pearson_r = if (sd(z$proportion.NNLS) > 0 && sd(z$proportion.ridge_NNLS) > 0)
               cor(z$proportion.NNLS, z$proportion.ridge_NNLS) else NA_real_,
             n_samples = nrow(z))
}))
write_tsv(agreement, file.path(out_dir, "deconvolution_algorithm_agreement.tsv"))

nnls_wide <- reshape(props[props$method == "NNLS", c("sample_id", "cell_type", "proportion")],
                     idvar = "sample_id", timevar = "cell_type", direction = "wide")
prop_mat <- as.matrix(nnls_wide[, grepl("^proportion", names(nnls_wide)), drop = FALSE])
variable_columns <- apply(prop_mat, 2, sd) > 1e-8
if (sum(variable_columns) >= 2) {
  comp_pc <- prcomp(prop_mat[, variable_columns, drop = FALSE], scale. = TRUE)$x[, 1:2, drop = FALSE]
} else if (sum(variable_columns) == 1) {
  comp_pc <- cbind(PC1 = as.numeric(scale(prop_mat[, variable_columns])), PC2 = 0)
  rownames(comp_pc) <- rownames(prop_mat)
} else {
  comp_pc <- matrix(0, nrow(prop_mat), 2,
                    dimnames = list(rownames(prop_mat), c("PC1", "PC2")))
}
model_data <- merge(axis, qc, by = c("sample_id", "patient_id", "pair_id", "cohort", "side"))
model_data$composition_PC1 <- comp_pc[match(model_data$sample_id, nnls_wide$sample_id), 1]
model_data$log_library_size <- log10(model_data$library_size)
model_data$side <- factor(model_data$side, levels = c("convex", "concave"))
extract_side <- function(data, formula, label) {
  fit <- lm(formula, data)
  co <- summary(fit)$coefficients
  data.frame(cohort = data$cohort[1], model = label,
             estimate_concave_minus_convex = co["sideconcave", "Estimate"],
             se = co["sideconcave", "Std. Error"],
             p_value = co["sideconcave", "Pr(>|t|)"], residual_df = df.residual(fit))
}
model_sensitivity <- do.call(rbind, lapply(split(model_data, model_data$cohort), function(z) rbind(
  extract_side(z, score ~ factor(pair_id) + side, "biology_only"),
  extract_side(z, score ~ factor(pair_id) + side + log_library_size, "biology_plus_library_QC"),
  extract_side(z, score ~ factor(pair_id) + side + composition_PC1, "biology_plus_composition_PC1")
)))
write_tsv(model_data, file.path(out_dir, "frozen_axis_model_data.tsv"))
write_tsv(model_sensitivity, file.path(out_dir, "frozen_axis_model_sensitivity.tsv"))

# Leave-one-patient-out within each diagnosis.
loo <- do.call(rbind, lapply(split(axis_eff, axis_eff$cohort), function(z) {
  do.call(rbind, lapply(z$pair_id, function(left_out) {
    kept <- z[z$pair_id != left_out, ]
    data.frame(cohort = z$cohort[1], omitted_pair = left_out,
               n_retained = nrow(kept), mean_paired_difference = mean(kept$paired_difference),
               positive_retained = sum(kept$paired_difference > 0))
  }))
}))
write_tsv(loo, file.path(out_dir, "frozen_axis_leave_one_patient_out.tsv"))

# Frozen pathway projections, exact sign-flip tests, and AIS-versus-CS contrasts.
frozen <- readRDS("data/processed/frozen_pathway_definitions.rds")
path_scores <- gene_set_zscore(logcpm_all, frozen$pathways, min_genes = 5)
path_long <- do.call(rbind, lapply(rownames(path_scores), function(gs) data.frame(
  gene_set = gs, sample_id = colnames(path_scores), score = path_scores[gs, ],
  patient_id = info_all$patient_id[match(colnames(path_scores), info_all$processed_sample_id)],
  pair_id = info_all$pair_id[match(colnames(path_scores), info_all$processed_sample_id)],
  cohort = info_all$cohort_label[match(colnames(path_scores), info_all$processed_sample_id)],
  side = as.character(info_all$side[match(colnames(path_scores), info_all$processed_sample_id)])
)))
path_eff <- merge(path_long[path_long$side == "concave", ], path_long[path_long$side == "convex", ],
                  by = c("gene_set", "pair_id", "patient_id", "cohort"),
                  suffixes = c(".concave", ".convex"))
path_eff$paired_difference <- path_eff$score.concave - path_eff$score.convex
path_tests <- do.call(rbind, lapply(split(path_eff, interaction(path_eff$cohort, path_eff$gene_set)), function(z) {
  cbind(cohort = z$cohort[1], gene_set = z$gene_set[1], exact_sign_flip(z$paired_difference))
}))
path_tests$FDR_exact_within_cohort <- ave(path_tests$exact_sign_flip_p, path_tests$cohort,
                                         FUN = function(p) p.adjust(p, "BH"))

path_interaction <- do.call(rbind, lapply(split(path_eff, path_eff$gene_set), function(z) {
  ais <- z$paired_difference[z$cohort == "AIS"]
  cs <- z$paired_difference[z$cohort == "CS"]
  obs <- mean(ais) - mean(cs)
  all_values <- c(ais, cs)
  null <- apply(ais_combinations, 2, function(idx) mean(all_values[idx]) - mean(all_values[-idx]))
  data.frame(gene_set = z$gene_set[1], mean_AIS = mean(ais), mean_CS = mean(cs),
             AIS_minus_CS = obs,
             exact_patient_label_permutation_p = mean(abs(null) >= abs(obs) - 1e-12))
}))
path_interaction$FDR_interaction <- p.adjust(path_interaction$exact_patient_label_permutation_p, "BH")
write_tsv(path_long, file.path(out_dir, "frozen_pathway_joint_scores.tsv"))
write_tsv(path_eff, file.path(out_dir, "frozen_pathway_joint_paired_effects.tsv"))
write_tsv(path_tests, file.path(out_dir, "frozen_pathway_within_cohort_tests.tsv"))
write_tsv(path_interaction, file.path(out_dir, "frozen_pathway_AIS_CS_interaction_tests.tsv"))

# Machine-readable interpretation gate.
cs_test <- axis_tests[axis_tests$cohort == "CS", ]
summary <- data.frame(
  claim = c("CS_frozen_axis_direction", "AIS_CS_axis_difference", "CS_gene_level_DE",
            "AIS_by_side_gene_interaction", "disease_specificity_conclusion"),
  result = c(
    sprintf("%d/%d positive; exact P=%.6g", cs_test$positive_pairs, cs_test$n_pairs, cs_test$exact_sign_flip_p),
    sprintf("AIS-CS mean difference=%.6g; exact P=%.6g", axis_interaction$difference_in_means,
            axis_interaction$exact_patient_label_permutation_p),
    sprintf("%d genes at FDR<0.05", sum(de_cs$FDR < 0.05, na.rm = TRUE)),
    sprintf("%d genes at FDR<0.05", sum(de_int$FDR < 0.05, na.rm = TRUE)),
    if (cs_test$exact_sign_flip_p < 0.05) "axis_reproduces_in_CS_not_AIS_specific"
    else if (axis_interaction$exact_patient_label_permutation_p < 0.05) "axis_differs_between_AIS_and_CS"
    else "inconclusive_no_evidence_of_specificity"
  ),
  role = c("negative_control", "formal_specificity_test", "negative_control",
           "formal_specificity_test", "interpretation_gate")
)
write_tsv(summary, file.path(out_dir, "specificity_summary.tsv"))
capture_versions(c("readxl", "edgeR", "limma", "Matrix", "Seurat", "nnls"),
                 "logs/gse254300_cs_specificity_R_versions.tsv")
message("GSE254300 CS disease-specificity control complete for five paired patients.")
