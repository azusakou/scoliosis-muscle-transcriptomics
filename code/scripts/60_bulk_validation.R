#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readxl); library(edgeR); library(limma); library(Matrix); library(Seurat); library(nnls)
})
source("scripts/lib_common.R")
set.seed(254300)

raw_root <- raw_data_root()
bulk_file <- file.path(raw_root, "bulk", "GSE254300_readcount.xlsx")
audit <- read.delim("metadata/audit/bulk_samples.tsv", check.names = FALSE)
audit <- audit[audit$cohort == "AIS", ]
validate_paired_samples(audit, "processed_sample_id", 5L)
raw <- as.data.frame(read_excel(bulk_file))
stopifnot(all(audit$processed_sample_id %in% names(raw)))
raw$GeneID <- make.unique(as.character(raw$GeneID))
counts <- as.matrix(raw[, audit$processed_sample_id, drop = FALSE])
rownames(counts) <- raw$GeneID
storage.mode(counts) <- "numeric"
info <- audit[match(colnames(counts), audit$processed_sample_id), ]
info$side <- factor(info$side, levels = c("convex", "concave"))

keep <- filterByExpr(counts, group = info$side, min.count = 10)
y <- DGEList(counts[keep, , drop = FALSE])
y <- calcNormFactors(y)
design <- model.matrix(~factor(patient_id) + side, info)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "sideconcave")
de <- topTags(qlf, n = Inf, sort.by = "none")$table
de$gene <- rownames(de); de$contrast <- "concave_minus_convex"; de$n_paired_patients <- 5L
write_tsv(de, "results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz")
logcpm <- cpm(y, log = TRUE, prior.count = 2)
write_tsv(data.frame(gene = rownames(logcpm), logcpm), "results/tables/bulk/GSE254300_AIS_TMM_logCPM.tsv.gz")
projection_matrix <- label_blind_logcpm(counts, y$samples$norm.factors)
projection_logcpm <- projection_matrix$logcpm
write_tsv(data.frame(gene = rownames(projection_logcpm), projection_logcpm),
          "results/tables/bulk/GSE254300_AIS_label_blind_projection_logCPM.tsv.gz")

qc <- data.frame(sample_id = colnames(counts), patient_id = info$patient_id, side = as.character(info$side),
                 library_size = colSums(counts), detected_genes = colSums(counts > 0),
                 normalization_factor = y$samples$norm.factors)
pc <- prcomp(t(logcpm), scale. = FALSE)
qc$PC1 <- pc$x[qc$sample_id, 1]; qc$PC2 <- pc$x[qc$sample_id, 2]
write_tsv(qc, "results/tables/bulk/GSE254300_AIS_QC.tsv")

# Project the exact frozen pathway definitions into bulk.
frozen <- readRDS("data/processed/frozen_pathway_definitions.rds")
path_scores <- gene_set_zscore(logcpm, frozen$pathways, min_genes = 5)
path_long <- do.call(rbind, lapply(rownames(path_scores), function(gs) data.frame(
  gene_set = gs, sample_id = colnames(path_scores), score = path_scores[gs, ],
  patient_id = info$patient_id[match(colnames(path_scores), info$processed_sample_id)],
  side = as.character(info$side[match(colnames(path_scores), info$processed_sample_id)]))))
path_effect <- merge(path_long[path_long$side == "concave", ], path_long[path_long$side == "convex", ],
                     by = c("gene_set", "patient_id"), suffixes = c(".concave", ".convex"))
path_effect$paired_difference <- path_effect$score.concave - path_effect$score.convex
path_test <- do.call(rbind, lapply(split(path_effect, path_effect$gene_set), function(z) {
  tt <- t.test(z$paired_difference)
  data.frame(gene_set = z$gene_set[1], mean_paired_difference = mean(z$paired_difference),
             se = sd(z$paired_difference)/sqrt(nrow(z)), p_value = tt$p.value, n_pairs = nrow(z))
}))
path_test$FDR <- p.adjust(path_test$p_value, method = "BH")
write_tsv(path_long, "results/tables/bulk/frozen_pathway_scores.tsv")
write_tsv(path_effect, "results/tables/bulk/frozen_pathway_paired_effects.tsv")
write_tsv(path_test, "results/tables/bulk/frozen_pathway_paired_tests.tsv")

# Project the discovery-only signature. The validation side labels and paired
# effects are not used to choose genes; only label-blind CPM detectability is
# applied to accommodate the bulk feature space.
signature <- read.delim("data/processed/frozen_multicellular_signature.tsv")
projection <- project_frozen_signature(projection_logcpm, signature)
axis_score <- projection$score
usage <- projection$weights
usage$dataset <- "GSE254300_AIS"
usage$detectability_rule <- sprintf("CPM>=%.1f_in_at_least_%d_libraries_label_blind",
                                    projection_matrix$min_cpm,
                                    projection_matrix$min_samples)
write_tsv(usage, "results/tables/bulk/frozen_multicellular_signature_projection_genes.tsv")
axis <- data.frame(sample_id = names(axis_score), score = as.numeric(axis_score),
                   patient_id = info$patient_id[match(names(axis_score), info$processed_sample_id)],
                   side = as.character(info$side[match(names(axis_score), info$processed_sample_id)]))
axis_eff <- merge(axis[axis$side == "concave", ], axis[axis$side == "convex", ], by = "patient_id", suffixes = c(".concave", ".convex"))
axis_eff$paired_difference <- axis_eff$score.concave - axis_eff$score.convex
write_tsv(axis, "results/tables/bulk/multicellular_axis_scores.tsv")
write_tsv(axis_eff, "results/tables/bulk/multicellular_axis_paired_effects.tsv")

# Cell composition sensitivity: a scRNA reference and two non-negative solvers.
obj <- readRDS("data/processed/scrna/scrna_annotated.rds")
obj <- JoinLayers(obj, assay = "SoupX")
sc <- LayerData(obj, assay = "SoupX", layer = "counts")
ct <- as.character(obj$cell_type)
eligible <- names(which(table(ct) >= 50))
reference <- sapply(eligible, function(k) Matrix::rowSums(sc[, ct == k, drop = FALSE]))
reference <- sweep(reference, 2, colSums(reference), "/") * 1e6
reference_log <- log2(reference + 1)
specificity <- do.call(rbind, lapply(seq_along(eligible), function(j) {
  other <- rowMeans(reference_log[, -j, drop = FALSE])
  data.frame(gene = rownames(reference_log), cell_type = eligible[j], specificity = reference_log[, j] - other,
             reference_logCPM = reference_log[, j])
}))
specificity <- specificity[!grepl("^(MT-|RPS|RPL|HBA|HBB)", specificity$gene), ]
markers <- do.call(rbind, lapply(split(specificity, specificity$cell_type), function(z) head(z[order(z$specificity, decreasing = TRUE), ], 40)))
markers <- markers[markers$gene %in% rownames(logcpm), ]
marker_genes <- unique(markers$gene)
ref_use <- reference_log[marker_genes, eligible, drop = FALSE]
bulk_use <- logcpm[marker_genes, , drop = FALSE]
deconv_one <- function(b, method) {
  if (method == "NNLS") x <- nnls::nnls(ref_use, b)$x else {
    lambda <- 1
    fn <- function(x) sum((b - ref_use %*% x)^2) + lambda * sum(x^2)
    x <- optim(rep(1/ncol(ref_use), ncol(ref_use)), fn, method = "L-BFGS-B", lower = 0)$par
  }
  if (sum(x) > 0) x / sum(x) else x
}
props <- do.call(rbind, lapply(c("NNLS", "ridge_NNLS"), function(method) do.call(rbind, lapply(colnames(bulk_use), function(s) {
  x <- deconv_one(bulk_use[, s], method)
  data.frame(method = method, sample_id = s, cell_type = eligible, proportion = x)
}))))
props$patient_id <- info$patient_id[match(props$sample_id, info$processed_sample_id)]
props$side <- as.character(info$side[match(props$sample_id, info$processed_sample_id)])
write_tsv(markers, "results/tables/bulk/deconvolution_marker_reference.tsv")
write_tsv(data.frame(
  reference_assay = "SoupX_ambient_corrected_counts",
  marker_selection = "top40_logCPM_specificity_per_eligible_cell_type",
  minimum_cells_per_cell_type = 50L,
  n_marker_genes = length(marker_genes)
), "results/tables/bulk/deconvolution_reference_manifest.tsv")
write_tsv(props, "results/tables/bulk/deconvolution_two_methods.tsv")
wide_props <- reshape(props, idvar = c("sample_id", "patient_id", "side", "cell_type"), timevar = "method", direction = "wide")
agreement <- do.call(rbind, lapply(split(wide_props, wide_props$cell_type), function(z) data.frame(
  cell_type = z$cell_type[1], pearson_r = if (sd(z$proportion.NNLS) > 0 && sd(z$proportion.ridge_NNLS) > 0)
    cor(z$proportion.NNLS, z$proportion.ridge_NNLS) else NA_real_, n_samples = nrow(z))))
write_tsv(agreement, "results/tables/bulk/deconvolution_algorithm_agreement.tsv")

# Biology-only, QC-only, biology+QC, and composition-adjusted score models.
model_data <- merge(axis, qc, by = c("sample_id", "patient_id", "side"))
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
  comp_pc <- matrix(0, nrow(prop_mat), 2, dimnames = list(rownames(prop_mat), c("PC1", "PC2")))
}
model_data$composition_PC1 <- comp_pc[match(model_data$sample_id, nnls_wide$sample_id), 1]
model_data$composition_PC2 <- comp_pc[match(model_data$sample_id, nnls_wide$sample_id), 2]
model_data$log_library_size <- log10(model_data$library_size)
extract_side <- function(formula, label) {
  fit <- lm(formula, model_data)
  co <- summary(fit)$coefficients
  if (!"sideconvex" %in% rownames(co)) return(data.frame(model = label, estimate_concave_minus_convex = NA, se = NA, p_value = NA, residual_df = df.residual(fit)))
  data.frame(model = label, estimate_concave_minus_convex = -co["sideconvex", "Estimate"], se = co["sideconvex", "Std. Error"],
             p_value = co["sideconvex", "Pr(>|t|)"], residual_df = df.residual(fit))
}
model_data$side <- factor(model_data$side, levels = c("concave", "convex"))
models <- rbind(
  extract_side(score ~ factor(patient_id) + side, "biology_only"),
  data.frame(model = "QC_only_no_side_term", estimate_concave_minus_convex = NA, se = NA, p_value = NA,
             residual_df = df.residual(lm(score ~ log_library_size + PC1 + PC2, model_data))),
  extract_side(score ~ factor(patient_id) + side + log_library_size, "biology_plus_library_QC"),
  extract_side(score ~ factor(patient_id) + side + composition_PC1, "biology_plus_composition_PC1")
)
write_tsv(model_data, "results/tables/bulk/multicellular_axis_model_data.tsv")
write_tsv(models, "results/tables/bulk/multicellular_axis_model_sensitivity.tsv")
capture_versions(c("readxl", "edgeR", "limma", "nnls", "Seurat"), "logs/bulk_package_versions.tsv")
message("Bulk paired validation complete for five AIS patient pairs.")
