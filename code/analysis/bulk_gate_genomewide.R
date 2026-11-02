#!/usr/bin/env Rscript
# Fixed retrospective RNA-gate sensitivity; never selects a new signature or gene universe.
# Portable inputs are public-source count TSVs plus archived original A coefficients.
# Rscript bulk_gate_genomewide.R --input-dir=PATH --output-dir=PATH
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly = TRUE)
allowed <- c("--input-dir", "--output-dir")
flags <- sub("=.*$", "", args)
if (any(!flags %in% allowed) || any(!grepl("=.+$", args)) || anyDuplicated(flags)) stop("Unknown, duplicate or empty CLI flag")
value <- function(flag, default) {v <- args[startsWith(args, paste0(flag, "="))]; if (length(v)) sub(paste0("^", flag, "="), "", v[1]) else default}
input <- value("--input-dir", "analysis/bulk_gate_inputs")
out <- value("--output-dir", "analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(digits = 17)
read_tab <- function(name) read.delim(if (grepl("gz$", name)) gzfile(file.path(input, name)) else file.path(input, name), check.names = FALSE)
read_matrix <- function(name) {d <- read_tab(name); m <- as.matrix(d[, -1, drop = FALSE]); rownames(m) <- d[[1]]; storage.mode(m) <- "double"; m}
write_tab <- function(x, suffix, gzip = FALSE) {
  name <- paste0("bulk_gate_", suffix, ".tsv", if (gzip) ".gz" else "")
  con <- if (gzip) gzfile(file.path(out, name), "wt") else file.path(out, name)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  if (inherits(con, "connection")) close(con)
}
sha <- function(p) digest::digest(file = p, algo = "sha256")
bindings <- read_tab("source_bindings.tsv")
stopifnot(!anyDuplicated(bindings$artifact))
current_hash <- vapply(bindings$artifact, function(n) sha(file.path(input, n)), character(1))
stopifnot(identical(unname(current_hash), as.character(bindings$SHA256)))
A <- read_matrix("GSE254300_tested_counts.tsv.gz"); B <- read_matrix("GSE311180_tested_counts.tsv.gz")
oldA <- read_tab("original_A_gene_results.tsv.gz"); oldB <- read_tab("original_B_gene_results.tsv.gz")
common <- read_tab("common_gene_universe.tsv")$gene
map <- read_tab("samples.tsv"); rna <- read_tab("RNA_expected_pairs.tsv")
keep <- sort(rna$patient_id[rna$interpretation == "consistent_RNA_profiles"])
omit <- sort(rna$patient_id[rna$interpretation != "consistent_RNA_profiles"])
stopifnot(identical(keep, c("STB-050", "STB-062", "STB-085", "STB-086", "STB-095")),
          identical(omit, c("STB-084", "STB-100")), length(common) == 12995L,
          identical(rownames(A), oldA$gene), identical(rownames(B), oldB$gene))
ma <- map[map$cohort == "GSE254300", ]; mb7 <- map[map$cohort == "GSE311180", ]
stopifnot(identical(colnames(A), ma$sample_id), identical(colnames(B), mb7$sample_id),
          all(table(ma$patient_id, ma$side) == 1), all(table(mb7$patient_id, mb7$side) == 1),
          all(is.finite(A)), all(A >= 0), all(is.finite(B)), all(B >= 0))
fit_model <- function(y, info) {
  info$side <- factor(info$side, levels = c("convex", "concave"))
  design <- model.matrix(~factor(patient_id) + side, info)
  yy <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(yy, design, robust = TRUE)
  q <- glmQLFTest(fit, coef = "sideconcave")
  tab <- topTags(q, n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  list(y = yy, fit = fit, result = tab, info = info)
}
qc <- read_tab("original_B_QC.tsv")
yA <- calcNormFactors(DGEList(A)); recalculatedB7 <- calcNormFactors(DGEList(B))
yB7 <- DGEList(B)
# Primary factors are the archived original factors, exactly as the fixed plan
# states. Recomputing factors from decimal gene sums can perturb TMM rank ties
# at binary last bits; that diagnostic is preserved rather than called exact.
yB7$samples$norm.factors <- qc$normalization_factor[match(colnames(B), qc$sample_id)]
baselineA <- fit_model(yA, ma); baselineB <- fit_model(yB7, mb7)
audit <- data.frame(check = c("original_A_logFC", "original_B_logFC_archived_TMM", "original_B_factor_transfer",
                             "recalculated_B7_factor_vs_archived_diagnostic"),
  maximum_absolute_error = c(max(abs(baselineA$result$logFC - oldA$logFC)),
    max(abs(baselineB$result$logFC - oldB$logFC)),
    max(abs(yB7$samples$norm.factors - qc$normalization_factor[match(colnames(B), qc$sample_id)])),
    max(abs(recalculatedB7$samples$norm.factors - yB7$samples$norm.factors))),
  tolerance = c(1e-7, 1e-7, 1e-14, NA), guard_applies = c(TRUE, TRUE, TRUE, FALSE))
# The published audit is an original-run receipt, not an error measurement that
# should be silently replaced with the current operating system's receipt.
if (!file.exists(file.path(out,"bulk_gate_input_reproduction_audit.tsv")))
  write_tab(audit, "input_reproduction_audit")
stopifnot(all(audit$maximum_absolute_error[audit$guard_applies] < audit$tolerance[audit$guard_applies]))
nullA <- read_matrix("A32_logFC.tsv.gz"); gridA <- as.matrix(read_tab("A32_grid.tsv")[, -1, drop = FALSE])
nullA <- nullA[common, , drop = FALSE]
obsA <- oldA$logFC[match(common, oldA$gene)]
stopifnot(ncol(nullA) == 32L, nrow(gridA) == 32L,
          max(abs(nullA[, which(rowSums(gridA == 1) == 5)] - obsA)) < 1e-7)
mb <- mb7[mb7$patient_id %in% keep, ]; idx <- match(mb$sample_id, mb7$sample_id)
locked <- yB7[, idx, keep.lib.sizes = TRUE]
native_retmm <- calcNormFactors(DGEList(B[, idx, drop = FALSE]))
# Published estimated factors are explicit model inputs. Native TMM near-tie
# trimming can differ across libm implementations; diagnose it rather than
# pretending that a fresh native estimation is bit-identical on every OS.
parameter <- read_tab("published_retained_TMM_factors.tsv")
stopifnot(!anyDuplicated(parameter$sample_id), setequal(parameter$sample_id, mb$sample_id))
retmm <- DGEList(B[, idx, drop = FALSE])
retmm$samples$norm.factors <- parameter$retained_library_TMM_factor[match(mb$sample_id, parameter$sample_id)]
stopifnot(all(is.finite(retmm$samples$norm.factors)), all(retmm$samples$norm.factors > 0))
diagdir <- "logs/gate_runtime_diagnostics"
dir.create(diagdir, recursive = TRUE, showWarnings = FALSE)
write.table(data.frame(sample_id=mb$sample_id,
  published_retained_TMM_factor=retmm$samples$norm.factors,
  native_reestimated_TMM_factor=native_retmm$samples$norm.factors,
  native_minus_published=native_retmm$samples$norm.factors-retmm$samples$norm.factors),
  file.path(diagdir,"native_retained_TMM_diagnostic.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(audit,file.path(diagdir,"current_input_reproduction_audit.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
capture.output(sessionInfo(),file=file.path(diagdir,"sessionInfo.log"))
write_tab(mb, "subset_metadata")
write_tab(data.frame(sample_id = mb$sample_id, patient_id = mb$patient_id, side = mb$side,
  locked_original_norm_factor = locked$samples$norm.factors, retained_library_TMM_factor = retmm$samples$norm.factors,
  retained_tested_library_size = locked$samples$lib.size), "normalization_factors")
write_tab(data.frame(original_tested_B_genes = nrow(B), fixed_common_genes = length(common),
  subset_all_zero_tested_genes = sum(rowSums(B[, idx, drop = FALSE]) == 0),
  subset_all_zero_common_genes = sum(rowSums(B[common, idx, drop = FALSE]) == 0)), "fixed_universe_availability")
rank_unit <- function(x) {r <- rank(x, ties.method = "average"); r <- r - mean(r); r / sqrt(sum(r^2))}
ra <- apply(nullA, 2, rank_unit); sa <- sign(nullA)
gridB <- as.matrix(expand.grid(rep(list(c(-1, 1)), 5)))
write_tab(data.frame(configuration = 1:32, gridB, check.names = FALSE), "B32_grid")
summary <- list(); candidates <- list(); loo <- list(); looA <- list(); patient_genes <- list()
nominees <- c("ANGPT1", "LEP")
# The interval inverts the same moderated QL test at candidate H0 values by
# shifting the side offset; full-model means/deviance and empirical-Bayes QL
# moderation stay fixed. P at H0=0 must reproduce the published zero-null test.
ql_interval <- function(model, gene) {
  full <- model$fit
  gene_idx <- which(rownames(full) == gene)
  cidx <- match("sideconcave", colnames(full$design)); side <- full$design[, cidx]
  mle <- full$unshrunk.coefficients[gene_idx, cidx] / log(2)
  test_at <- function(delta) {
    # Preserve ALL genes: glmQLFTest caps denominator degrees of freedom by
    # sum(df.residual) across the complete fit. A one-gene slice would silently
    # change that cap and therefore change the original moderated test.
    shifted <- full
    shifted$offset <- as.matrix(full$offset)
    shifted$offset[gene_idx, ] <- shifted$offset[gene_idx, ] + side * delta * log(2)
    shifted$coefficients[gene_idx, cidx] <- full$coefficients[gene_idx, cidx] - delta * log(2)
    shifted$unshrunk.coefficients[gene_idx, cidx] <- full$unshrunk.coefficients[gene_idx, cidx] - delta * log(2)
    as.numeric(glmQLFTest(shifted, coef = cidx)$table$PValue[gene_idx])
  }
  expected <- model$result$PValue[match(gene, model$result$gene)]
  stopifnot(abs(test_at(0) - expected) < 1e-8, test_at(mle) > .95)
  step <- 1
  while (test_at(mle - step) > .05 && step < 1024) step <- step * 2
  lower <- uniroot(function(d) test_at(d) - .05, c(mle - step, mle), tol = 1e-7)$root
  step <- 1
  while (test_at(mle + step) > .05 && step < 1024) step <- step * 2
  upper <- uniroot(function(d) test_at(d) - .05, c(mle, mle + step), tol = 1e-7)$root
  stopifnot(lower <= mle, upper >= mle, abs(test_at(lower) - .05) < 1e-6,
            abs(test_at(upper) - .05) < 1e-6)
  c(MLE_unshrunk_logFC = mle, QL95_lower = lower, QL95_upper = upper,
    P_reproduction_absolute_error = abs(test_at(0) - expected))
}
for (label in c("locked_original_TMM", "retained_library_reTMM")) {
  y <- if (label == "locked_original_TMM") locked else retmm
  observed <- fit_model(y, mb)
  write_tab(observed$result, paste0("observed_genes_", label), TRUE)
  effect <- observed$result$logFC[match(common, observed$result$gene)]
  observed_rho <- cor(obsA, effect, method = "spearman")
  observed_sign <- mean(sign(obsA) == sign(effect))
  nullB <- matrix(NA_real_, length(common), 32, dimnames = list(common, paste0("configuration_", sprintf("%03d", 1:32))))
  for (j in which(gridB[, 1] == 1)) {
    sign <- gridB[j, ]; info <- mb
    flip <- sign[match(info$patient_id, keep)] < 0
    info$side[flip] <- ifelse(info$side[flip] == "concave", "convex", "concave")
    perm <- fit_model(y, info)
    nullB[, j] <- perm$result$logFC[match(common, perm$result$gene)]
    complement <- which(apply(gridB, 1, function(z) all(z == -sign)))
    nullB[, complement] <- -nullB[, j]
    cat(label, "B configuration", j, "complete\n")
  }
  stopifnot(max(abs(nullB[, 32] - effect)) < 1e-7)
  write_tab(data.frame(gene = common, nullB, check.names = FALSE), paste0("B32_logFC_", label), TRUE)
  rb <- apply(nullB, 2, rank_unit); sb <- sign(nullB)
  rho <- t(ra) %*% rb
  agree <- ((t(sa) %*% sb + t(abs(sa)) %*% abs(sb)) / 2 + t(sa == 0) %*% (sb == 0)) / length(common)
  p <- c(mean(abs(rho) >= abs(observed_rho) - 1e-12), mean(abs(agree - .5) >= abs(observed_sign - .5) - 1e-12))
  summary[[label]] <- data.frame(normalization = label, endpoint = c("spearman_logFC", "sign_concordance"),
    observed = c(observed_rho, observed_sign), exact_two_sided_p = p, BH_within_two_endpoint_family = p.adjust(p, "BH"),
    n_A_pairs = 5, n_B_pairs = 5, n_common_genes = length(common), n_joint_side_configurations = 1024)
  write_tab(data.frame(A_configuration = rep(1:32, 32), B_configuration = rep(1:32, each = 32),
    spearman_rho = as.numeric(rho), sign_concordance = as.numeric(agree)), paste0("exact_null_", label))
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  for (gene in nominees) {
    row <- observed$result[observed$result$gene == gene, ]
    stopifnot(nrow(row) == 1)
    d <- vapply(keep, function(pid) unname(logcpm[gene, mb$sample_id[mb$patient_id == pid & mb$side == "concave"]] -
      logcpm[gene, mb$sample_id[mb$patient_id == pid & mb$side == "convex"]]), numeric(1))
    ci <- ql_interval(observed, gene)
    tt <- t.test(d)
    candidates[[paste(label, gene)]] <- data.frame(normalization = label, gene = gene, fitted_shrunk_logFC = row$logFC,
      as.list(ci), QL_F = row$F, QL_P = row$PValue, genome_wide_FDR_original_17961 = row$FDR,
      patient_paired_logCPM_mean = mean(d), patient_paired_logCPM_SD = sd(d),
      patient_t_CI95_lower = tt$conf.int[1], patient_t_CI95_upper = tt$conf.int[2],
      positive_pairs = sum(d > 0), n_pairs = 5)
    patient_genes[[paste(label, gene)]] <- data.frame(normalization = label, gene = gene,
      patient_id = keep, paired_logCPM_difference = d)
  }
  if (label == "locked_original_TMM") for (pid in keep) {
    retained <- mb$patient_id != pid
    model <- fit_model(y[, retained, keep.lib.sizes = TRUE], mb[retained, ])
    e <- model$result$logFC[match(common, model$result$gene)]
    loo[[pid]] <- data.frame(omitted_B_patient = pid, n_B_pairs = 4, n_common_genes = length(common),
      spearman_rho = cor(obsA, e, method = "spearman"), sign_concordance = mean(sign(obsA) == sign(e)),
      role = "descriptive_patient_influence_not_independent_replication")
  }
  if (label == "locked_original_TMM") for (pid in sort(unique(ma$patient_id))) {
    retained <- ma$patient_id != pid
    model <- fit_model(yA[, retained, keep.lib.sizes = TRUE], ma[retained, ])
    e <- model$result$logFC[match(common, model$result$gene)]
    looA[[pid]] <- data.frame(omitted_A_patient = pid, n_A_pairs = 4, n_B_pairs = 5,
      n_common_genes = length(common), spearman_rho = cor(e, effect, method = "spearman"),
      sign_concordance = mean(sign(e) == sign(effect)),
      role = "later_amended_descriptive_patient_influence_not_independent_replication")
  }
}
all_summary <- do.call(rbind, summary)
all_summary$BH_full_four_endpoint_sensitivity_family <- p.adjust(all_summary$exact_two_sided_p, "BH")
write_tab(all_summary, "genomewide_exact_summary")
candidate_table <- do.call(rbind, candidates)
candidate_table$BH_two_fixed_candidates_within_normalization <- ave(candidate_table$QL_P, candidate_table$normalization,
  FUN = function(p) p.adjust(p, "BH"))
candidate_table$BH_four_candidate_normalization_sensitivity_family <- p.adjust(candidate_table$QL_P, "BH")
write_tab(candidate_table, "ANGPT1_LEP_fixed_candidate_results")
write_tab(do.call(rbind, patient_genes), "fixed_candidate_patient_effects")
write_tab(do.call(rbind, loo), "B_patient_omission_influence")
write_tab(do.call(rbind, looA), "A_patient_omission_influence")
write_tab(data.frame(package = c("R", "edgeR"), version = c(as.character(getRversion()),
  as.character(packageVersion("edgeR")))), "runtime_versions")
write_tab(bindings, "verified_input_bindings")
cat("Fixed RNA-gate genome-wide and two-candidate sensitivity completed.\n")
