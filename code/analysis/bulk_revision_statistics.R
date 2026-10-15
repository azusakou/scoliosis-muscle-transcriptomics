#!/usr/bin/env Rscript
source("scripts/analysis_paths.R")
# Reviewer-requested analyses. Run from the repository root:
# Rscript analysis/bulk_revision_statistics.R
# This script reads original outputs and writes only analysis/bulk_* files.
suppressPackageStartupMessages(library(metafor))
out <- "analysis"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_tab <- function(x, name) write.table(x, file.path(out, paste0("bulk_", name, ".tsv")),
                                        sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
read_matrix <- function(path) {
  x <- read.delim(gzfile(path), check.names = FALSE)
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "double"
  rownames(m) <- x[[1]]; m
}
sign_grid <- function(n) as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
stats_one <- function(x, label) {
  x <- as.numeric(x); n <- length(x); se <- sd(x) / sqrt(n)
  null <- as.numeric(sign_grid(n) %*% x / n)
  wt <- suppressWarnings(wilcox.test(x, mu = 0, alternative = "two.sided", exact = TRUE,
                                    correct = FALSE))
  ci <- mean(x) + c(-1, 1) * qt(.975, n - 1) * se
  data.frame(analysis = label, n_pairs = n, mean_paired_difference = mean(x),
             paired_difference_SD = sd(x), SE = se, CI95_lower = ci[1], CI95_upper = ci[2],
             cohen_dz = mean(x) / sd(x), median_paired_difference = median(x),
             positive_pairs = sum(x > 0), t_test_p = t.test(x)$p.value,
             exact_sign_flip_p = mean(abs(null) >= abs(mean(x)) - 1e-12),
             wilcoxon_V = unname(wt$statistic), exact_wilcoxon_p = wt$p.value,
             n_sign_configurations = length(null), min_attainable_sign_flip_p = 2 / 2^n)
}
pooled <- read.delim("results/tables/meta_validation_requant/frozen_axis_12_patient_effects.tsv")
cohort_summaries <- do.call(rbind, lapply(split(pooled, pooled$cohort), function(z)
  stats_one(z$paired_difference, z$cohort[1])))
summaries <- rbind(cohort_summaries, stats_one(pooled$paired_difference, "pooled_12_AIS"))
write_tab(summaries, "patient_effect_statistics")
loo <- do.call(rbind, lapply(seq_len(nrow(pooled)), function(i) {
  cbind(omitted_cohort = pooled$cohort[i], omitted_patient = pooled$patient_id[i],
        stats_one(pooled$paired_difference[-i], "pooled_leave_one_patient_out"))
}))
write_tab(loo, "pooled_leave_one_patient_out")
write_tab(data.frame(configuration = seq_len(2^nrow(pooled)),
                     permuted_mean = as.numeric(sign_grid(nrow(pooled)) %*%
                                                 pooled$paired_difference / nrow(pooled))),
          "pooled_exact_sign_flip_distribution")

# Both patient-weighted and equal-cohort estimands are explicit. A cohort factor
# permits different cohort means; this model assumes a common within-cohort variance.
pooled$cohort <- factor(pooled$cohort)
fit <- lm(paired_difference ~ 0 + cohort, pooled)
coefs <- coef(fit); vc <- vcov(fit)
cohort_weights <- as.numeric(table(pooled$cohort)) / nrow(pooled)
model_contrast <- function(w, label) {
  est <- sum(w * coefs); se <- sqrt(as.numeric(t(w) %*% vc %*% w)); df <- df.residual(fit)
  ci <- est + c(-1, 1) * qt(.975, df) * se
  data.frame(analysis = label, estimate = est, SE = se, CI95_lower = ci[1],
             CI95_upper = ci[2], residual_df = df,
             t_test_p = 2 * pt(-abs(est / se), df),
             cohort_weight_GSE254300 = w[1], cohort_weight_GSE311180 = w[2])
}
write_tab(rbind(model_contrast(cohort_weights, "cohort_factor_patient_weighted"),
                model_contrast(c(.5, .5), "cohort_factor_equal_cohort_weighted")),
          "cohort_stratified_estimates")
hetero_test <- t.test(paired_difference ~ cohort, pooled)
write_tab(data.frame(contrast = "GSE254300_minus_GSE311180", estimate = diff(-hetero_test$estimate),
                     CI95_lower = hetero_test$conf.int[1], CI95_upper = hetero_test$conf.int[2],
                     welch_df = unname(hetero_test$parameter), p_value = hetero_test$p.value),
          "cohort_mean_difference")
cs <- cohort_summaries
meta <- rma.uni(yi = cs$mean_paired_difference, vi = cs$SE^2,
                method = "REML", test = "knha", slab = cs$analysis)
fixed <- rma.uni(yi = cs$mean_paired_difference, vi = cs$SE^2, method = "FE")
meta_row <- function(m, label) data.frame(analysis = label, n_cohorts = m$k,
  estimate = as.numeric(m$b), SE = m$se, CI95_lower = m$ci.lb, CI95_upper = m$ci.ub,
  p_value = m$pval, Q = m$QE, Q_df = m$k - 1, heterogeneity_p = m$QEp,
  I2_percent = m$I2, tau2 = m$tau2,
  limitation = "two_cohorts_cannot_precisely_estimate_between_study_heterogeneity")
write_tab(rbind(meta_row(meta, "REML_Hartung_Knapp"), meta_row(fixed, "inverse_variance_fixed_effect")),
          "meta_analysis_heterogeneity")

# The stored paired scores remain locked in the first exclusion analysis.
# A second analysis re-z-scores the retained cohort libraries while holding the
# label-blind detectable universe and frozen gene identities/signs unchanged.
matA <- read_matrix("results/tables/bulk/GSE254300_AIS_label_blind_projection_logCPM.tsv.gz")
matB <- read_matrix("results/tables/bulk_gse311180_requant/GSE311180_label_blind_projection_logCPM.tsv.gz")
sig <- read.delim("results/tables/bulk/frozen_multicellular_signature.tsv")
weights <- aggregate(weight ~ gene, sig, mean)
weights <- weights[is.finite(weights$weight) & weights$weight != 0, ]
score <- function(m, w = weights) {
  w <- w[w$gene %in% rownames(m), ]
  z <- t(scale(t(m[w$gene, , drop = FALSE]))); z[!is.finite(z)] <- 0
  colSums(z * w$weight) / sum(abs(w$weight))
}
metaA <- read.delim("metadata/audit/bulk_samples.tsv")
metaA <- metaA[metaA$cohort == "AIS", ]; metaA$sample_id <- metaA$processed_sample_id
metaB <- read.delim("metadata/audit/gse311180_samples.tsv", check.names = FALSE)
metaB <- metaB[metaB$sample_id %in% colnames(matB), ]
paired <- function(s, meta) {
  vapply(sort(unique(meta$patient_id)), function(id)
    unname(s[meta$sample_id[meta$patient_id == id & meta$side == "concave"]] -
           s[meta$sample_id[meta$patient_id == id & meta$side == "convex"]]), numeric(1))
}
effA <- paired(score(matA), metaA); effB <- paired(score(matB), metaB)
checkA <- pooled[pooled$cohort == "GSE254300", ]
checkB <- pooled[pooled$cohort == "GSE311180_requant", ]
stopifnot(max(abs(effA[checkA$patient_id] - checkA$paired_difference)) < 1e-12,
          max(abs(effB[checkB$patient_id] - checkB$paired_difference)) < 1e-12)
age20 <- unique(metaB$patient_id[metaB$age == 20])
stopifnot(identical(age20, "STB-050"))
retainedB <- metaB[!metaB$patient_id %in% age20, ]
newB <- paired(score(matB[, retainedB$sample_id, drop = FALSE]), retainedB)
age_results <- rbind(
  stats_one(pooled$paired_difference[!pooled$patient_id %in% age20], "exclude_age20_locked_scores_pooled11"),
  stats_one(effB[!names(effB) %in% age20], "exclude_age20_locked_scores_GSE311180_6"),
  stats_one(c(effA, newB), "exclude_age20_rez_retained_libraries_pooled11"),
  stats_one(newB, "exclude_age20_rez_retained_libraries_GSE311180_6"))
age_results$omitted_patient <- age20
age_results$site_caveat <- "STB-050_concave_T6-T7_convex_T8-T9_age_and_site_exclusion_are_same_sensitivity"
write_tab(age_results, "age20_and_site_exclusion")

# Post hoc consequence of the fixed RNA pairing gates; no relabelling or retuning.
rna_pairs <- read.delim(analysis_input("analysis/audit_gse311180_rna_expected_pairs.tsv"))
omitRNA <- sort(rna_pairs$patient_id[rna_pairs$interpretation != "consistent_RNA_profiles"])
stopifnot(length(omitRNA)==2, all(omitRNA %in% names(effB)))
retainedRNA <- metaB[!metaB$patient_id %in% omitRNA, ]
rezRNA <- paired(score(matB[, retainedRNA$sample_id, drop=FALSE]), retainedRNA)
lockedRNA <- effB[!names(effB) %in% omitRNA]
rna_results <- rbind(
 stats_one(c(effA, lockedRNA), "RNA_gate_exclusion_locked_pooled10"),
 stats_one(lockedRNA, "RNA_gate_exclusion_locked_GSE311180_5"),
 stats_one(c(effA, rezRNA), "RNA_gate_exclusion_rez_pooled10"),
 stats_one(rezRNA, "RNA_gate_exclusion_rez_GSE311180_5"))
rna_results$omitted_patients <- paste(omitRNA, collapse=";")
rna_results$scope <- "post_hoc_fixed_signature_sensitivity_not_anatomical_validation_or_new_primary_endpoint"
write_tab(rna_results, "RNA_pairing_gate_exclusion")
write_tab(rbind(data.frame(patient_id=names(lockedRNA),scaling="locked",paired_difference=lockedRNA),
 data.frame(patient_id=names(rezRNA),scaling="retained_libraries_rez",paired_difference=rezRNA)),
 "RNA_pairing_gate_exclusion_patient_effects")

# Joint expression-decile matched random signature: each random gene and frozen
# signed slot is shared by both cohorts. Joint bins also preserve availability
# patterns and therefore the exact 122/126 gene score denominators. No target
# gene is eligible for a random set. One set has unique gene identities.
bins <- function(m) {
  x <- rank(rowMeans(m), ties.method = "average")
  b <- cut(x, breaks = quantile(x, probs = seq(0, 1, .1)),
           include.lowest = TRUE, labels = FALSE); names(b) <- rownames(m); b
}
binA <- bins(matA); binB <- bins(matB)
universe <- union(names(binA), names(binB))
key <- paste(ifelse(is.na(binA[universe]), 0, binA[universe]),
             ifelse(is.na(binB[universe]), 0, binB[universe]), sep = "_")
names(key) <- universe
slots <- weights[weights$gene %in% universe, ]
slots$joint_bin <- key[slots$gene]
pools <- split(universe[!universe %in% weights$gene], key[!universe %in% weights$gene])
bin_manifest <- do.call(rbind, lapply(unique(slots$joint_bin), function(k)
  data.frame(joint_bin = k, target_genes = sum(slots$joint_bin == k),
             candidate_genes = length(pools[[k]]))))
stopifnot(all(bin_manifest$candidate_genes >= bin_manifest$target_genes))
write_tab(bin_manifest, "pooled_random_signature_joint_bin_manifest")
gene_mean_z_effect <- function(m, meta) {
  z <- t(scale(t(m))); z[!is.finite(z)] <- 0
  pairs <- sort(unique(meta$patient_id))
  d <- vapply(pairs, function(id) z[, meta$sample_id[meta$patient_id == id & meta$side == "concave"]] -
                z[, meta$sample_id[meta$patient_id == id & meta$side == "convex"]], numeric(nrow(z)))
  setNames(rowMeans(d), rownames(m))
}
gA <- gene_mean_z_effect(matA, metaA); gB <- gene_mean_z_effect(matB, metaB)
useA <- slots$gene %in% names(gA); useB <- slots$gene %in% names(gB)
denA <- sum(abs(slots$weight[useA])); denB <- sum(abs(slots$weight[useB]))
set.seed(20261012)
B <- 10000L
null <- matrix(NA_real_, B, 3, dimnames = list(NULL, c("GSE254300_mean", "GSE311180_mean", "pooled12_mean")))
first_set <- NULL
draw_gene_sets <- character(B)
for (i in seq_len(B)) {
  random_genes <- character(nrow(slots))
  for (k in unique(slots$joint_bin)) {
    idx <- which(slots$joint_bin == k)
    random_genes[idx] <- sample(pools[[k]], length(idx), replace = FALSE)
  }
  stopifnot(!anyDuplicated(random_genes))
  draw_gene_sets[i] <- paste(sort(random_genes), collapse = "|")
  null[i, 1] <- sum(gA[random_genes[useA]] * slots$weight[useA]) / denA
  null[i, 2] <- sum(gB[random_genes[useB]] * slots$weight[useB]) / denB
  null[i, 3] <- (5 * null[i, 1] + 7 * null[i, 2]) / 12
  if (i == 1) first_set <- data.frame(target_gene = slots$gene, random_gene = random_genes,
                                     frozen_signed_slot = slots$weight, joint_bin = slots$joint_bin)
}
stopifnot(!anyDuplicated(draw_gene_sets))
write_tab(data.frame(n_random_draws = B, n_distinct_gene_sets = length(unique(draw_gene_sets)),
                     all_genes_unique_within_each_set = TRUE, seed = 20261012),
          "pooled_random_signature_uniqueness_audit")
observed <- c(mean(effA), mean(effB), mean(c(effA, effB)))
random_summary <- data.frame(analysis = colnames(null), n_random_signatures = B,
  observed_mean = observed, null_mean = colMeans(null), null_SD = apply(null, 2, sd),
  empirical_two_sided_p = (1 + colSums(abs(null) >= rep(abs(observed), each = B))) / (B + 1),
  n_available_nonzero_genes_GSE254300 = sum(useA), n_available_nonzero_genes_GSE311180 = sum(useB),
  n_shared_available_genes = sum(useA & useB), seed = 20261012,
  null_recipe = "same_unique_random_genes_and_frozen_signs_in_both_cohorts_joint_expression_deciles_and_detectability_patterns")
write_tab(random_summary, "pooled_matched_random_signature_summary")
write_tab(cbind(iteration = seq_len(B), as.data.frame(null)), "pooled_matched_random_signature_null")
write_tab(first_set, "pooled_random_signature_first_draw_audit")

# Diagnosis-by-side contrast is between patient-level paired effects. Welch CI
# estimates the interaction without assuming equal variances. A planning MDE is
# derived from the observed pooled SD of paired differences; it is not post-hoc
# power or evidence of equivalence. Normal-model exact permutation power is also
# simulated, with all 252 balanced patient assignments used within each trial.
joint <- read.delim("results/tables/disease_control_gse254300_cs/frozen_axis_joint_paired_effects.tsv")
a <- joint$paired_difference[joint$cohort == "AIS"]
c <- joint$paired_difference[joint$cohort == "CS"]
tt <- t.test(a, c, var.equal = FALSE)
sigma <- sqrt(((length(a) - 1) * var(a) + (length(c) - 1) * var(c)) / (length(a) + length(c) - 2))
write_tab(rbind(stats_one(a, "joint_normalized_AIS5"), stats_one(c, "joint_normalized_CS5")),
          "interaction_within_diagnosis_statistics")
write_tab(data.frame(contrast = "AIS_minus_CS_patient_paired_effect", difference = mean(a) - mean(c),
                     SE_Welch = sqrt(var(a) / 5 + var(c) / 5), CI95_lower = tt$conf.int[1],
                     CI95_upper = tt$conf.int[2], welch_df = unname(tt$parameter), welch_p = tt$p.value,
                     exact_patient_label_permutation_p = read.delim("results/tables/disease_control_gse254300_cs/frozen_axis_AIS_CS_interaction_test.tsv")$exact_patient_label_permutation_p,
                     SD_AIS = sd(a), SD_CS = sd(c), pooled_SD_paired_differences = sigma,
                     n_AIS_pairs = 5, n_CS_pairs = 5, n_label_assignments = choose(10, 5),
                     minimum_two_sided_exact_p = 2 / choose(10, 5)), "interaction_uncertainty")
parametric_mde <- power.t.test(n = 5, sd = sigma, power = .8, sig.level = .05,
                              type = "two.sample", alternative = "two.sided", strict = TRUE)$delta
grid_d <- seq(0, 3, .1)
set.seed(20261013)
nsim <- 20000L
base <- matrix(rnorm(10 * nsim), 10, nsim)
assignments <- combn(10, 5)
contrast <- matrix(-1 / 5, choose(10, 5), 10)
for (i in seq_len(nrow(contrast))) contrast[i, assignments[, i]] <- 1 / 5
perm_power <- function(d) {
  x <- base; x[1:5, ] <- x[1:5, ] + d
  observed_d <- colMeans(x[1:5, , drop = FALSE]) - colMeans(x[6:10, , drop = FALSE])
  null_d <- contrast %*% x
  p <- colMeans(abs(null_d) >= rep(abs(observed_d) - 1e-12, each = nrow(null_d)))
  mean(p <= .05)
}
sim_power <- vapply(grid_d, perm_power, numeric(1))
power_table <- data.frame(standardized_interaction_d = grid_d, interaction_score_units = grid_d * sigma,
                         exact_permutation_power_estimate = sim_power,
                         Monte_Carlo_SE = sqrt(sim_power * (1 - sim_power) / nsim), simulations = nsim,
                         parametric_t_power = vapply(grid_d, function(d) power.t.test(n = 5,
                           delta = d * sigma, sd = sigma, sig.level = .05, type = "two.sample",
                           alternative = "two.sided", strict = TRUE)$power, numeric(1)))
write_tab(power_table, "interaction_power_curve")
sim_mde_d <- approx(sim_power, grid_d, xout = .8, ties = mean)$y
write_tab(data.frame(alpha = .05, target_power = .8, n_patients_per_diagnosis = 5,
  planning_SD = sigma, parametric_t_MDE_score_units = parametric_mde,
  parametric_t_MDE_standardized = parametric_mde / sigma,
  exact_permutation_MDE_score_units = sim_mde_d * sigma,
  exact_permutation_MDE_standardized = sim_mde_d, simulations_per_grid_point = nsim,
  simulation_seed = 20261013,
  assumptions = "independent_normally_distributed_patient_paired_differences_equal_variance_fixed_SD_estimate_no_multiple_testing_adjustment",
  interpretation = "planning_sensitivity_not_posthoc_power_and_not_equivalence_test"), "interaction_minimum_detectable_effect")

# Patient-level exact side permutations of TMM logCPM summaries retain every
# gene's co-expression and every patient's paired dependence. They do not use
# genes as biological replicates and are labelled separately from edgeR logFC.
tmA <- read_matrix("results/tables/bulk/GSE254300_AIS_TMM_logCPM.tsv.gz")
tmB <- read_matrix("results/tables/bulk_gse311180_requant/GSE311180_TMM_logCPM.tsv.gz")
replication <- read.delim(gzfile("results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz"))
common <- replication$gene
stopifnot(all(common %in% rownames(tmA)), all(common %in% rownames(tmB)))
pair_matrix <- function(m, meta) {
  vapply(sort(unique(meta$patient_id)), function(id)
    m[common, meta$sample_id[meta$patient_id == id & meta$side == "concave"]] -
      m[common, meta$sample_id[meta$patient_id == id & meta$side == "convex"]], numeric(length(common)))
}
dA <- pair_matrix(tmA, metaA); dB <- pair_matrix(tmB, metaB)
sgA <- sign_grid(5); sgB <- sign_grid(7)
meanA <- dA %*% t(sgA) / 5; meanB <- dB %*% t(sgB) / 7
rank_unit <- function(x) { r <- rank(x, ties.method = "average"); r <- r - mean(r); r / sqrt(sum(r^2)) }
rA <- apply(meanA, 2, rank_unit); rB <- apply(meanB, 2, rank_unit)
rho <- t(rA) %*% rB
sA <- sign(meanA); sB <- sign(meanB)
agree <- ((t(sA) %*% sB + t(abs(sA)) %*% abs(sB)) / 2 +
            t(sA == 0) %*% (sB == 0)) / length(common)
obsrho <- cor(rowMeans(dA), rowMeans(dB), method = "spearman")
obsagree <- mean(sign(rowMeans(dA)) == sign(rowMeans(dB)))
write_tab(data.frame(analysis = "patient_side_exact_TMM_logCPM_mean_effects", n_common_genes = length(common),
                    observed_spearman_rho = obsrho, observed_sign_concordance = obsagree,
                    exact_two_sided_rho_p = mean(abs(rho) >= abs(obsrho) - 1e-12),
                    exact_two_sided_sign_concordance_p = mean(abs(agree - .5) >= abs(obsagree - .5) - 1e-12),
                    n_joint_side_configurations = 4096,
                    original_edgeR_spearman_rho = cor(replication$logFC.GSE254300,
                                                      replication$logFC.GSE311180_requant, method = "spearman"),
                    original_edgeR_sign_concordance = mean(replication$sign_concordant),
                    null_limitation = "logCPM_effect_statistic_not_refitted_edgeR_logFC"), "gene_concordance_patient_permutation_summary")
write_tab(data.frame(GSE254300_configuration = rep(seq_len(32), 128),
                    GSE311180_configuration = rep(seq_len(128), each = 32),
                    spearman_rho = as.numeric(rho), sign_concordance = as.numeric(agree)),
          "gene_concordance_patient_permutation_null")
gmeanA <- rowMeans(dA); gmeanB <- rowMeans(dB)
gvarA <- apply(dA, 1, var) / 5; gvarB <- apply(dB, 1, var) / 7
Q <- (gmeanA - gmeanB)^2 / (gvarA + gvarB)
qh <- data.frame(gene = common, mean_paired_logCPM_GSE254300 = gmeanA,
                 mean_paired_logCPM_GSE311180 = gmeanB, sampling_var_GSE254300 = gvarA,
                 sampling_var_GSE311180 = gvarB, Q = Q, Q_df = 1,
                 Q_p = pchisq(Q, 1, lower.tail = FALSE), I2_percent = pmax(0, (Q - 1) / Q) * 100)
qh$heterogeneity_FDR <- p.adjust(qh$Q_p, "BH")
write_tab(qh, "gene_heterogeneity_logCPM_descriptive")
siggenes <- weights$gene
toprep <- replication[replication$FDR.GSE311180_requant < .05, ]
write_tab(data.frame(common_tested_genes = nrow(replication),
  GSE311180_FDR05_common_genes = nrow(toprep), same_sign_in_GSE254300 = sum(toprep$sign_concordant),
  frozen_nonzero_signature_overlap = sum(toprep$gene %in% siggenes),
  non_signature_genes = sum(!toprep$gene %in% siggenes),
  non_signature_same_sign = sum(toprep$sign_concordant[!toprep$gene %in% siggenes]),
  nominal_Q_p05_genes = sum(qh$Q_p < .05), heterogeneity_FDR05_genes = sum(qh$heterogeneity_FDR < .05)),
  "gene_replication_signature_overlap")
write_tab(data.frame(package = c("R", "metafor"), version = c(as.character(getRversion()),
  as.character(packageVersion("metafor")))), "statistics_runtime_versions")
cat("Reviewer-requested bulk statistics written to analysis/bulk_*.tsv\n")
