#!/usr/bin/env Rscript
# Reproduce archived random-signature recipes and remove the two zero-weight
# conflicts from the GSE254300 random-signature target size. Both original
# observed scores are verified from saved label-blind projection matrices.
out <- "analysis"
write_tab <- function(x, suffix) write.table(x, file.path(out, paste0("bulk_", suffix, ".tsv")),
                                          sep = "\t", quote = FALSE, row.names = FALSE)
rd <- function(path) {a <- read.delim(gzfile(path), check.names = FALSE); m <- as.matrix(a[, -1]); rownames(m) <- a[[1]]; m}
sig <- read.delim("results/tables/bulk/frozen_multicellular_signature.tsv")
allw <- aggregate(weight ~ gene, sig, mean)
A <- rd("results/tables/bulk/GSE254300_AIS_label_blind_projection_logCPM.tsv.gz")
B <- rd("results/tables/bulk_gse311180_requant/GSE311180_label_blind_projection_logCPM.tsv.gz")
ma <- read.delim("metadata/audit/bulk_samples.tsv"); ma <- ma[ma$cohort == "AIS", ]; ma$sample_id <- ma$processed_sample_id
mb <- read.delim("metadata/audit/gse311180_samples.tsv"); mb <- mb[mb$sample_id %in% colnames(B), ]
gene_effect <- function(m, meta) {
  z <- t(scale(t(m))); z[!is.finite(z)] <- 0
  x <- vapply(sort(unique(meta$patient_id)), function(pid) z[, meta$sample_id[meta$patient_id == pid & meta$side == "concave"]] -
                z[, meta$sample_id[meta$patient_id == pid & meta$side == "convex"]], numeric(nrow(z)))
  setNames(rowMeans(x), rownames(m))
}
ga <- gene_effect(A, ma); gb <- gene_effect(B, mb)
meanA <- rowMeans(A); meanB <- rowMeans(B)
ba <- cut(meanA, quantile(meanA, seq(0, 1, .1)), include.lowest = TRUE, labels = FALSE); names(ba) <- rownames(A)
rb <- rank(meanB, ties.method = "average")
bb <- cut(rb, quantile(rb, seq(0, 1, .1)), include.lowest = TRUE, labels = FALSE); names(bb) <- rownames(B)
draw <- function(w, bins) {
  target_bins <- bins[w$gene]; x <- character(nrow(w))
  for (bin in unique(target_bins)) {
    idx <- which(target_bins == bin)
    pool <- names(bins)[bins == bin & !names(bins) %in% w$gene]
    stopifnot(length(pool) >= length(idx))
    x[idx] <- sample(pool, length(idx), replace = FALSE)
  }
  stopifnot(!anyDuplicated(x)); x
}
run <- function(g, bins, w, n, seed, shuffle, label) {
  set.seed(seed)
  observed <- sum(g[w$gene] * w$weight) / sum(abs(w$weight))
  null <- vapply(seq_len(n), function(i) {
    gs <- draw(w, bins)
    sw <- if (shuffle) sample(w$weight) else w$weight
    sum(g[gs] * sw) / sum(abs(sw))
  }, numeric(1))
  list(summary = data.frame(analysis = label, n_sampled_gene_slots = nrow(w),
    n_nonzero_signed_slots = sum(w$weight != 0), n_random_sets = n, seed = seed,
    signs_shuffled_globally_per_set = shuffle, observed_mean = observed,
    empirical_two_sided_p = (1 + sum(abs(null) >= abs(observed))) / (n + 1),
    null_mean = mean(null), null_SD = sd(null)), null = data.frame(iteration = seq_len(n), random_mean = null))
}
wAold <- allw[allw$gene %in% names(ga), ]; wA <- wAold[wAold$weight != 0, ]
wB <- allw[allw$gene %in% names(gb) & allw$weight != 0, ]
oldA <- run(ga, ba, wAold, 1000, 700722, FALSE, "archived_GSE254300_recipe_124_slots_122_nonzero")
oldB <- run(gb, bb, wB, 1000, 1368770, TRUE, "archived_GSE311180_recipe_126_nonzero_slots")
stopifnot(abs(oldA$summary$empirical_two_sided_p - read.delim("results/tables/robustness/bulk_axis_matched_random_gene_sets_summary.tsv")$empirical_two_sided_p) < 1e-12,
          abs(oldB$summary$empirical_two_sided_p - read.delim("results/tables/bulk_gse311180_requant/frozen_axis_random_signature_test.tsv")$expression_decile_and_size_matched_empirical_p) < 1e-12)
newA <- run(ga, ba, wA, 10000, 700722, FALSE, "GSE254300_corrected_122_nonzero_slots_original_bin_and_sign_recipe")
newB <- run(gb, bb, wB, 10000, 1368770, TRUE, "GSE311180_126_nonzero_slots_original_bin_and_sign_recipe_10k")
write_tab(rbind(oldA$summary, oldB$summary, newA$summary, newB$summary), "cohort_random_null_recipe_audit")
write_tab(newA$null, "GSE254300_corrected_random_signature_null")
write_tab(newB$null, "GSE311180_expanded_random_signature_null")
cat("Archived cohort nulls reproduced; corrected 122/126 slot nulls written.\n")
