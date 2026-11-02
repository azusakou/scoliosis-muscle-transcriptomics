#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(edgeR); library(limma) })
source("scripts/lib_common.R")
set.seed(700722)

axis <- read.delim("results/tables/bulk/multicellular_axis_paired_effects.tsv")
obs <- mean(axis$paired_difference)
signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(axis))))
perm_means <- as.numeric(signs %*% axis$paired_difference / nrow(axis))
exact_p <- mean(abs(perm_means) >= abs(obs) - 1e-12)
write_tsv(data.frame(observed_mean = obs, exact_two_sided_p = exact_p, n_label_configurations = nrow(signs)),
          "results/tables/robustness/bulk_axis_exact_side_label_permutation.tsv")
write_tsv(data.frame(permutation_id = seq_along(perm_means), permuted_mean = perm_means),
          "results/tables/robustness/bulk_axis_permutation_distribution.tsv")

loo <- do.call(rbind, lapply(axis$patient_id, function(pid) data.frame(
  omitted_patient = pid, mean_paired_difference = mean(axis$paired_difference[axis$patient_id != pid]),
  sign_preserved = sign(mean(axis$paired_difference[axis$patient_id != pid])) == sign(obs))))
write_tsv(loo, "results/tables/robustness/bulk_axis_leave_one_patient_out.tsv")

logcpm <- read.delim(gzfile("results/tables/bulk/GSE254300_AIS_label_blind_projection_logCPM.tsv.gz"), check.names = FALSE)
rownames(logcpm) <- logcpm$gene; logcpm$gene <- NULL; logcpm <- as.matrix(logcpm)
meta <- read.delim("metadata/audit/bulk_samples.tsv")
meta <- meta[meta$cohort == "AIS", ]; meta <- meta[match(colnames(logcpm), meta$processed_sample_id), ]
sig <- read.delim("data/processed/frozen_multicellular_signature.tsv")
weights <- aggregate(weight ~ gene, sig, mean)
weights <- weights[weights$gene %in% rownames(logcpm), ]
gene_mean <- rowMeans(logcpm)
bins <- cut(gene_mean, quantile(gene_mean, probs = seq(0, 1, .1), na.rm = TRUE), include.lowest = TRUE, labels = FALSE)
names(bins) <- rownames(logcpm)
score_gene_set <- function(genes, w) {
  z <- t(scale(t(logcpm[genes, , drop = FALSE]))); z[!is.finite(z)] <- 0
  s <- colSums(z * w) / sum(abs(w))
  mean(vapply(unique(meta$patient_id), function(pid) {
    s[meta$patient_id == pid & meta$side == "concave"] - s[meta$patient_id == pid & meta$side == "convex"]
  }, numeric(1)))
}
# Archived 1,000-draw, 124-slot null includes two zero-weight conflicts.
# Superseded for inference by analysis/bulk_cohort_matched_random_null.R.
random_effects <- replicate(1000, {
  sampled <- sample_matched_genes_without_replacement(
    weights$gene, bins[weights$gene], bins, weights$gene)
  score_gene_set(sampled, weights$weight)
})
random_p <- (1 + sum(abs(random_effects) >= abs(obs))) / (1 + length(random_effects))
write_tsv(data.frame(observed = obs, empirical_two_sided_p = random_p, n_random_sets = length(random_effects),
                     matching = "signature_size_and_bulk_mean_expression_decile"),
          "results/tables/robustness/bulk_axis_matched_random_gene_sets_summary.tsv")
write_tsv(data.frame(iteration = seq_along(random_effects), random_effect = random_effects,
                     n_unique_genes = nrow(weights), target_unique_genes = nrow(weights)),
          "results/tables/robustness/bulk_axis_matched_random_gene_sets.tsv.gz")

weights <- weights[weights$weight != 0, , drop = FALSE]
# Count the actual nonzero scoring genes in this sensitivity summary.
high_cut <- quantile(gene_mean, .99)
keep_sig <- weights$gene[gene_mean[weights$gene] < high_cut]
high_removed <- score_gene_set(keep_sig, weights$weight[match(keep_sig, weights$gene)])
write_tsv(data.frame(analysis = c("all_signature_genes", "remove_top_1pct_mean_expression"),
                     mean_paired_difference = c(obs, high_removed), n_genes = c(nrow(weights), length(keep_sig))),
          "results/tables/robustness/bulk_axis_high_expression_sensitivity.tsv")

# scRNA pseudobulk leave-one-donor effects are the two donor-specific estimates themselves.
sc <- read.delim("results/tables/programs/scrna_pathway_paired_effects_by_donor.tsv")
sc_loo <- sc[, c("cell_type", "gene_set", "patient_id", "paired_difference")]
names(sc_loo)[3:4] <- c("retained_patient", "leave_one_donor_out_effect")
sc_loo$inferential_status <- "single_donor_sensitivity_only"
write_tsv(sc_loo, "results/tables/robustness/scrna_pathway_leave_one_donor_out.tsv")

message("Robustness and negative controls complete.")
