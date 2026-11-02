#!/usr/bin/env Rscript
# Refit the published full-cohort edgeR side configurations from retained counts.
# This recreates models; it does not realign reads or certify anatomical labels.
suppressPackageStartupMessages({library(edgeR); library(digest); library(jsonlite)})
if (as.character(getRversion()) != "4.4.2" || as.character(packageVersion("edgeR")) != "4.4.2")
  stop("Full retained-model replay requires the recorded R4.4.2/edgeR4.4.2 runtime")
root <- normalizePath(".", mustWork = TRUE)
input <- "analysis/bulk_gate_inputs"
bindings <- read.delim(file.path(input, "source_bindings.tsv"), check.names = FALSE)
stopifnot(!anyDuplicated(bindings$artifact),
          identical(unname(vapply(bindings$artifact, function(n)
            digest(file = file.path(input, n), algo = "sha256"), character(1))), as.character(bindings$SHA256)))
read_tab <- function(path) read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path, check.names = FALSE)
read_matrix <- function(path) {
  tab <- read_tab(path); value <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(value) <- "double"; rownames(value) <- tab[[1]]; value
}
common <- read_tab(file.path(input, "common_gene_universe.tsv"))$gene
samples <- read_tab(file.path(input, "samples.tsv"))
stopifnot(length(common) == 12995L, !anyDuplicated(common), !anyDuplicated(samples$sample_id))
dir.create("logs/full_concordance_refit", recursive = TRUE, showWarnings = FALSE)
run <- tempfile("run-", tmpdir = "logs/full_concordance_refit")
stopifnot(!file.exists(run)); dir.create(run)
rev <- file.path(run, "analysis"); dir.create(rev)
dir.create(file.path(rev, "bulk_gate_inputs"))
dir.create(file.path(run, "results/tables/meta_validation_requant"), recursive = TRUE)
replication_path <- "results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz"
stopifnot(file.copy(replication_path, file.path(run, replication_path), overwrite = FALSE))
replication <- read_tab(replication_path)
stopifnot(identical(common, replication$gene))
write_tsv <- function(tab, path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con)); write.table(tab, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_coefficients <- function(value, path) {
  # Preserve 17 significant digits, as in the original coefficient exporter.
  columns <- as.data.frame(lapply(seq_len(ncol(value)), function(i) sprintf("%.17g", value[, i])),
                           stringsAsFactors = FALSE)
  names(columns) <- sprintf("configuration_%03d", seq_len(ncol(value)))
  write_tsv(cbind(data.frame(gene = rownames(value)), columns), path)
}
fit_table <- function(y, info) {
  design <- model.matrix(~factor(patient_id) + side, info)
  fitted <- glmQLFit(estimateDisp(y, design, robust = TRUE), design, robust = TRUE)
  topTags(glmQLFTest(fitted, coef = "sideconcave"), n = Inf, sort.by = "none")$table
}
model_close <- function(actual, expected, label) {
  stopifnot(identical(dim(actual), dim(expected)), all(is.finite(actual)), all(is.finite(expected)))
  relative_error <- max(abs(actual - expected) / pmax(1, abs(expected)))
  if (relative_error > 1e-7) stop("Existing model-refit tolerance exceeded: ", label, " (", relative_error, ")")
  relative_error
}
fits <- list(); audits <- list()
for (cohort in c("GSE254300", "GSE311180")) {
  tag <- if (cohort == "GSE254300") "A" else "B"
  n <- if (tag == "A") 5L else 7L
  counts <- read_matrix(file.path(input, paste0(cohort, "_tested_counts.tsv.gz")))
  original <- read_tab(file.path(input, paste0("original_", tag, "_gene_results.tsv.gz")))
  info <- samples[samples$cohort == cohort, , drop = FALSE]
  stopifnot(identical(rownames(counts), original$gene), identical(colnames(counts), info$sample_id),
            length(unique(info$patient_id)) == n, all(table(info$patient_id, info$side) == 1),
            all(is.finite(counts)), all(counts >= 0))
  info$side <- factor(info$side, levels = c("convex", "concave"))
  y <- calcNormFactors(DGEList(counts))
  if (tag == "B") {
    qc <- read_tab(file.path(input, "original_B_QC.tsv"))
    stopifnot(setequal(qc$sample_id, info$sample_id), !anyDuplicated(qc$sample_id))
    y$samples$norm.factors <- qc$normalization_factor[match(info$sample_id, qc$sample_id)]
  }
  observed <- fit_table(y, info)
  expected <- as.matrix(original[, c("logFC", "logCPM", "F", "PValue", "FDR")])
  actual <- as.matrix(observed[original$gene, colnames(expected)])
  observed_error <- model_close(actual, expected, paste(cohort, "observed model"))
  stopifnot(identical(unname(sign(actual[, "logFC"])), unname(sign(expected[, "logFC"]))),
            identical(unname(actual[, "PValue"] < .05), unname(expected[, "PValue"] < .05)),
            identical(unname(actual[, "FDR"] < .05), unname(expected[, "FDR"] < .05)))
  patients <- sort(unique(info$patient_id))
  grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  coefficients <- matrix(NA_real_, length(common), nrow(grid), dimnames = list(common, NULL))
  for (j in which(grid[, 1] == 1)) {
    permuted <- info
    flipped <- grid[j, match(info$patient_id, patients)] < 0
    side <- as.character(info$side)
    side[flipped] <- ifelse(side[flipped] == "concave", "convex", "concave")
    permuted$side <- factor(side, levels = c("convex", "concave"))
    fitted <- fit_table(y, permuted)
    coefficients[, j] <- fitted[common, "logFC"]
    complement <- which(apply(grid, 1, function(z) all(z == -grid[j, ])))
    coefficients[, complement] <- -coefficients[, j]
    if (j %% 16 == 0 || j == max(which(grid[, 1] == 1)))
      cat(cohort, "refitted configuration", j, "of", nrow(grid), "\n")
  }
  reference_path <- if (tag == "A") file.path(input, "A32_logFC.tsv.gz") else
    "analysis/bulk_edgeR_null_logFC_GSE311180.tsv.gz"
  reference <- read_matrix(reference_path)
  stopifnot(identical(rownames(reference), common))
  colnames(coefficients) <- colnames(reference)
  error <- model_close(coefficients, reference, paste(cohort, "all coefficients"))
  stopifnot(identical(sign(coefficients), sign(reference)))
  coefficient_path <- if (tag == "A") file.path(rev, "bulk_gate_inputs/A32_logFC.tsv.gz") else
    file.path(rev, "bulk_edgeR_null_logFC_GSE311180.tsv.gz")
  write_coefficients(coefficients, coefficient_path)
  grid_table <- data.frame(configuration = seq_len(nrow(grid)), grid)
  names(grid_table)[-1] <- paste0("patient_sign_", patients)
  grid_path <- if (tag == "A") file.path(rev, "bulk_gate_inputs/A32_grid.tsv") else
    file.path(rev, "bulk_edgeR_null_grid_GSE311180.tsv")
  write_tsv(grid_table, grid_path)
  audits[[cohort]] <- data.frame(cohort = cohort, tested_genes = nrow(original),
    common_genes = length(common), configurations = nrow(grid), distinct_model_refits = nrow(grid)/2,
    max_scaled_observed_error = observed_error, max_scaled_coefficient_error = error)
  fits[[cohort]] <- coefficients
}
rank_unit <- function(x) {r <- rank(x, ties.method = "average"); r <- r - mean(r); r / sqrt(sum(r^2))}
rankA <- apply(fits$GSE254300, 2, rank_unit); rankB <- apply(fits$GSE311180, 2, rank_unit)
rho <- t(rankA) %*% rankB
sA <- sign(fits$GSE254300); sB <- sign(fits$GSE311180)
agree <- ((t(sA) %*% sB + t(abs(sA)) %*% abs(sB)) / 2 + t(sA == 0) %*% (sB == 0)) / length(common)
observed_rho <- rho[32, 128]; observed_sign <- agree[32, 128]
null <- data.frame(GSE254300_configuration = rep(1:32, 128),
                   GSE311180_configuration = rep(1:128, each = 32),
                   spearman_rho = as.numeric(rho), sign_concordance = as.numeric(agree))
summary <- data.frame(observed_spearman_rho = observed_rho, observed_sign_concordance = observed_sign,
  exact_two_sided_rho_p = mean(abs(rho) >= abs(observed_rho) - 1e-12),
  exact_two_sided_sign_concordance_p = mean(abs(agree - .5) >= abs(observed_sign - .5) - 1e-12),
  n_common_genes = length(common), n_joint_side_configurations = 4096)
reference_null <- read_tab("analysis/bulk_gene_concordance_edgeR_permutation_null.tsv")
reference_summary <- read_tab("analysis/bulk_gene_concordance_edgeR_permutation_summary.tsv")
stopifnot(identical(as.matrix(null[, 1:2]), as.matrix(reference_null[, 1:2])),
          max(abs(null$spearman_rho - reference_null$spearman_rho)) < 2e-10,
          max(abs(null$sign_concordance - reference_null$sign_concordance)) < 2e-10)
for (field in names(summary))
  stopifnot(abs(summary[[field]] - reference_summary[[field]]) < 2e-10)
write_tsv(null, file.path(rev, "bulk_gene_concordance_edgeR_permutation_null.tsv"))
write_tsv(summary, file.path(rev, "bulk_gene_concordance_edgeR_permutation_summary.tsv"))
write_tsv(do.call(rbind, audits), file.path(run, "model_refit_checks.tsv"))
capture.output(sessionInfo(), file = file.path(run, "sessionInfo.log"))
report <- list(status = "MODEL_REFITS_PASSED_PENDING_INDEPENDENT_CHECK", run_root = normalizePath(run),
  source_script_SHA256 = digest(file = "scripts/refit_full_concordance.R", algo = "sha256"),
  source_bindings_SHA256 = digest(file = file.path(input, "source_bindings.tsv"), algo = "sha256"),
  scientific_model_tolerance = 1e-7, concordance_arithmetic_tolerance = 2e-10,
  coefficient_signs_and_P_FDR_decisions_preserved = TRUE,
  refitted_configurations = c(32L, 128L), distinct_permutation_model_fits = 80L,
  joint_configurations = 4096L, normalization = "original_tested_universes_and_archived_B_TMM_factors",
  scope = "Retained public aggregate count refits; no new raw-read alignment or identity claim")
write_json(report, file.path(run, "model_refit_receipt.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(normalizePath(run), "logs/full_concordance_refit/latest_run.txt")
cat("Full original-cohort refit completed:", normalizePath(run), "\n")
