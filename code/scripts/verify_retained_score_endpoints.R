#!/usr/bin/env Rscript
# Independent downstream score tests from retained sample scores/expression.
# No GMT scoring, count-level CS interaction model or cell algorithm is rerun.
source("scripts/lib_common.R")
out <- "logs/retained_score_endpoints"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(path) read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path,
                              check.names = FALSE, stringsAsFactors = FALSE)
key <- function(x, columns) do.call(paste, c(x[columns], sep = "\034"))
checks <- list()
compare <- function(actual, path, keys, name) {
  expected <- rd(path)
  if (length(keys)) {
    a <- key(actual, keys); b <- key(expected, keys)
    stopifnot(!anyDuplicated(a), !anyDuplicated(b), setequal(a, b))
    actual <- actual[match(b, a), , drop = FALSE]
  }
  stopifnot(nrow(actual) == nrow(expected), all(names(expected) %in% names(actual)))
  for (column in names(expected)) {
    a <- actual[[column]]; b <- expected[[column]]
    if (!identical(is.na(a), is.na(b))) stop(name, ": missingness differs in ", column)
    keep <- !is.na(a)
    if (is.numeric(a) && is.numeric(b)) {
      if (any(!is.finite(a[keep])) || any(!is.finite(b[keep])) ||
          any(abs(a[keep] - b[keep]) > pmax(1e-12, 1e-10 * abs(b[keep]))))
        stop(name, ": arithmetic mismatch in ", column)
      if (grepl("^n_|^positive_pairs$", column) && any(a[keep] != b[keep]))
        stop(name, ": count mismatch in ", column)
      if (grepl("p_value|_p$|FDR|BH", column)) {
        if (any(a[keep] < 0 | a[keep] > 1 | b[keep] < 0 | b[keep] > 1))
          stop(name, ": probability outside [0,1] in ", column)
        if (any((a[keep] < .05) != (b[keep] < .05)))
          stop(name, ": changed probability decision")
      }
      if (column %in% c("mean_paired_difference", "median_paired_difference",
                        "AIS_minus_CS", "difference_in_means") &&
          any(sign(a[keep]) != sign(b[keep])))
        stop(name, ": changed effect sign in ", column)
    } else if (!identical(as.character(a[keep]), as.character(b[keep]))) {
      stop(name, ": label mismatch in ", column)
    }
  }
  checks[[length(checks) + 1L]] <<- data.frame(check = name, rows = nrow(expected),
    fields = paste(names(expected), collapse = ";"), status = "PASS", reference = path)
  write_tsv(actual[, names(expected), drop = FALSE], file.path(out, paste0(name, ".tsv")))
  invisible(actual)
}
paired <- function(data, keys, positive = "concave", negative = "convex", side = "side") {
  select <- c(keys, "score")
  p <- data[data[[side]] == positive, select, drop = FALSE]
  n <- data[data[[side]] == negative, select, drop = FALSE]
  stopifnot(!anyDuplicated(key(p, keys)), !anyDuplicated(key(n, keys)),
            setequal(key(p, keys), key(n, keys)))
  result <- merge(p, n, by = keys, suffixes = c(".concave", ".convex"))
  result$paired_difference <- result$score.concave - result$score.convex
  result
}
test <- function(x) {
  stopifnot(all(is.finite(x)), length(x) >= 2L)
  n <- length(x); mean <- mean(x); se <- sd(x) / sqrt(n)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  null <- as.numeric(signs %*% x / n)
  list(n_pairs = n, positive_pairs = sum(x > 0), mean_paired_difference = mean,
    median_paired_difference = median(x), se = se,
    p = 2 * pt(-abs(mean / se), n - 1),
    exact = mean(abs(null) >= abs(mean) - 1e-12), min_p = 2 / 2^n)
}
pathway_stats <- function(scores, style, keys = c("gene_set", "patient_id")) {
  effects <- paired(scores, keys)
  split_columns <- if ("cohort" %in% keys) c("cohort", "gene_set") else "gene_set"
  groups <- split(effects, key(effects, split_columns))
  result <- do.call(rbind, lapply(groups, function(z) {
    t <- test(z$paired_difference)
    label <- z[1, split_columns, drop = FALSE]
    if (style == "first") cbind(label, data.frame(mean_paired_difference = t$mean_paired_difference,
      se = t$se, p_value = t$p, n_pairs = t$n_pairs)) else if (style == "second")
      cbind(label, data.frame(n_pairs = t$n_pairs, mean_paired_difference = t$mean_paired_difference,
        positive_pairs = t$positive_pairs, t_test_p = t$p, exact_sign_flip_p = t$exact)) else
      cbind(label, data.frame(n_pairs = t$n_pairs, positive_pairs = t$positive_pairs,
        mean_paired_difference = t$mean_paired_difference, median_paired_difference = t$median_paired_difference,
        paired_t_p = t$p, exact_sign_flip_p = t$exact, min_attainable_two_sided_p = t$min_p))
  }))
  if (style == "first") result$FDR <- p.adjust(result$p_value, "BH")
  if (style == "second") {
    result$FDR_t_test <- p.adjust(result$t_test_p, "BH")
    result$FDR_exact <- p.adjust(result$exact_sign_flip_p, "BH")
  }
  if (style == "joint") result$FDR_exact_within_cohort <-
    ave(result$exact_sign_flip_p, result$cohort, FUN = function(p) p.adjust(p, "BH"))
  list(effects = effects, tests = result)
}

# Each complete retained pathway family is tested; no nominally favorable subset.
A <- pathway_stats(rd("results/tables/bulk/frozen_pathway_scores.tsv"), "first")
compare(A$tests, "results/tables/bulk/frozen_pathway_paired_tests.tsv", "gene_set", "first_cohort_pathway_tests")
B <- pathway_stats(rd("results/tables/bulk_gse311180_requant/frozen_pathway_scores.tsv"), "second")
compare(B$tests, "results/tables/bulk_gse311180_requant/frozen_pathway_paired_tests.tsv", "gene_set", "second_cohort_pathway_tests")
joint_scores <- rd("results/tables/disease_control_gse254300_cs/frozen_pathway_joint_scores.tsv")
J <- pathway_stats(joint_scores, "joint", c("gene_set", "pair_id", "patient_id", "cohort"))
compare(J$tests, "results/tables/disease_control_gse254300_cs/frozen_pathway_within_cohort_tests.tsv",
  c("cohort", "gene_set"), "joint_diagnosis_pathway_tests")
combinations <- combn(10L, 5L)
interaction_test <- function(ais, cs) {
  stopifnot(length(ais) == 5L, length(cs) == 5L)
  values <- c(ais, cs); observed <- mean(ais) - mean(cs)
  null <- apply(combinations, 2, function(i) mean(values[i]) - mean(values[-i]))
  list(observed = observed, null = null, p = mean(abs(null) >= abs(observed) - 1e-12))
}
interaction <- do.call(rbind, lapply(split(J$effects, J$effects$gene_set), function(z) {
  ais <- z$paired_difference[z$cohort == "AIS"]; cs <- z$paired_difference[z$cohort == "CS"]
  t <- interaction_test(ais, cs)
  data.frame(gene_set = z$gene_set[1], mean_AIS = mean(ais), mean_CS = mean(cs),
    AIS_minus_CS = t$observed, exact_patient_label_permutation_p = t$p)
}))
interaction$FDR_interaction <- p.adjust(interaction$exact_patient_label_permutation_p, "BH")
compare(interaction, "results/tables/disease_control_gse254300_cs/frozen_pathway_AIS_CS_interaction_tests.tsv",
  "gene_set", "diagnosis_pathway_interaction_tests")

# Contextual systems keep their original roles; they are not pooled as muscle pairs.
osteo <- pathway_stats(rd("results/tables/adjacent_tissue_gse324067/frozen_pathway_scores.tsv"), "second")
compare(osteo$tests, "results/tables/adjacent_tissue_gse324067/frozen_pathway_paired_tests.tsv",
  "gene_set", "osteoblast_context_pathway_tests")

# Rebuild the diagnosis-label test directly from the ten retained sample pairs.
axis <- rd("results/tables/disease_control_gse254300_cs/frozen_axis_joint_scores.tsv")
effects <- paired(axis, c("pair_id", "patient_id", "cohort"))
ais <- effects$paired_difference[effects$cohort == "AIS"]
cs <- effects$paired_difference[effects$cohort == "CS"]
t <- interaction_test(ais, cs)
summary <- data.frame(contrast = "AIS_minus_CS_paired_axis_difference", n_AIS_pairs = 5L, n_CS_pairs = 5L,
  mean_AIS = mean(ais), mean_CS = mean(cs), difference_in_means = t$observed,
  exact_patient_label_permutation_p = t$p, n_label_permutations = length(t$null))
compare(summary, "results/tables/disease_control_gse254300_cs/frozen_axis_AIS_CS_interaction_test.tsv",
  "contrast", "diagnosis_score_label_test")
compare(data.frame(permuted_difference = t$null),
  "results/tables/disease_control_gse254300_cs/frozen_axis_AIS_CS_interaction_null.tsv",
  character(), "diagnosis_score_all_252_labels")
within <- do.call(rbind, lapply(split(effects, effects$cohort), function(z) {
  t <- test(z$paired_difference)
  data.frame(cohort = z$cohort[1], n_pairs = t$n_pairs, positive_pairs = t$positive_pairs,
    mean_paired_difference = t$mean_paired_difference, median_paired_difference = t$median_paired_difference,
    paired_t_p = t$p, exact_sign_flip_p = t$exact, min_attainable_two_sided_p = t$min_p)
}))
compare(within, "results/tables/disease_control_gse254300_cs/frozen_axis_within_cohort_tests.tsv",
  "cohort", "joint_diagnosis_score_tests")

# The original second-cohort label-blind projection matrix and fixed weights
# suffice to recompute the exact original high-expression sensitivity recipe.
matrix <- rd("results/tables/bulk_gse311180_requant/GSE311180_label_blind_projection_logCPM.tsv.gz")
m <- as.matrix(matrix[, -1]); rownames(m) <- matrix$gene
signature <- rd("results/tables/bulk/frozen_multicellular_signature.tsv")
weights <- aggregate(weight ~ gene, signature, mean)
weights <- weights[weights$weight != 0 & weights$gene %in% rownames(m), ]
info <- rd("metadata/audit/gse311180_samples.tsv")
info <- info[match(colnames(m), info$sample_id), ]
stopifnot(!anyNA(info$sample_id), length(unique(info$patient_id)) == 7L, nrow(weights) == 126L)
project <- function(w) {
  z <- t(scale(t(m[w$gene, , drop = FALSE]))); z[!is.finite(z)] <- 0
  score <- colSums(z * w$weight) / sum(abs(w$weight))
  data.frame(sample_id = names(score), patient_id = info$patient_id, side = info$side, score = as.numeric(score))
}
original <- project(weights)
reference_scores <- rd("results/tables/bulk_gse311180_requant/frozen_axis_scores.tsv")
compare(original, "results/tables/bulk_gse311180_requant/frozen_axis_scores.tsv", "sample_id", "second_cohort_score_projection")
mean_expression <- rowMeans(m)
cutoff <- quantile(mean_expression, .99)
trimmed <- weights[mean_expression[weights$gene] < cutoff, ]
trimmed_scores <- project(trimmed)
trimmed_effects <- paired(trimmed_scores, "patient_id")
t <- test(trimmed_effects$paired_difference)
trimmed_test <- data.frame(n_signature_genes = nrow(weights), n_retained_genes = nrow(trimmed),
  mean_paired_difference = t$mean_paired_difference, positive_pairs = t$positive_pairs, t_test_p = t$p)
compare(trimmed_test, "results/tables/bulk_gse311180_requant/frozen_axis_remove_top1pct.tsv",
  character(), "second_cohort_original_top1pct_recipe")
write_tsv(trimmed_scores, file.path(out, "second_cohort_top1pct_sample_scores.tsv"))
write_tsv(trimmed_effects, file.path(out, "second_cohort_top1pct_patient_effects.tsv"))
write_tsv(trimmed, file.path(out, "second_cohort_top1pct_weights.tsv"))
write_tsv(do.call(rbind, checks), file.path(out, "checks.tsv"))
writeLines(c("Complete retained pathway-test families, diagnosis-label tests and original seven-pair top1% sensitivity recomputed.",
  "Reference results unchanged; P/FDR decisions and integer counts remain exact; arithmetic uses ABS1e-12/REL1e-10.",
  "Pathway sample scores and diagnosis sample scores remain declared inputs; no GMT scoring or CS count-level DE reconstruction."),
  file.path(out, "SCOPE.txt"))
cat("Retained score endpoint checks PASS:", length(checks), "complete table comparisons. See", out, "\n")
