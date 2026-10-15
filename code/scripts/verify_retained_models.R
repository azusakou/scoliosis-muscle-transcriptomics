#!/usr/bin/env Rscript
# Read-only checks from retained aggregate counts and model data.
# This does not reconstruct cell-level QC, ambient correction or tissue identity.
suppressPackageStartupMessages(library(edgeR))
close <- function(actual, expected, label, tolerance = 1e-7) {
  stopifnot(length(actual) == length(expected), identical(unname(is.na(actual)), unname(is.na(expected))))
  finite <- !is.na(actual)
  if (any(!is.finite(actual[finite])) || any(!is.finite(expected[finite])) ||
      any(abs(actual[finite] - expected[finite]) > tolerance * pmax(1, abs(expected[finite])))) {
    stop("Mismatch: ", label)
  }
}
original <- read.delim(gzfile("results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz"), check.names = FALSE)
checked_genes <- 0L
for (lineage in unique(original$cell_type)) {
  name <- gsub("[^A-Za-z0-9]+", "_", lineage)
  table <- read.delim(gzfile(paste0("results/tables/pseudobulk/counts/", name, ".tsv.gz")), check.names = FALSE)
  counts <- as.matrix(table[, -1]); rownames(counts) <- table$gene
  info <- data.frame(patient_id = sub("__.*", "", colnames(counts)),
                     side = factor(sub(".*__", "", colnames(counts)), levels = c("convex", "concave")))
  design <- model.matrix(~factor(patient_id) + side, info)
  keep <- filterByExpr(counts, design = design, min.count = 5)
  y <- calcNormFactors(DGEList(counts[keep, , drop = FALSE]))
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  actual <- topTags(glmQLFTest(fit, coef = "sideconcave"), n = Inf, sort.by = "none")$table
  expected <- original[original$cell_type == lineage, ]; rownames(expected) <- expected$gene
  stopifnot(setequal(rownames(actual), expected$gene))
  for (column in c("logFC", "logCPM", "F", "PValue", "FDR")) {
    close(actual[expected$gene, column], expected[[column]], paste(lineage, column))
  }
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  for (patient in unique(info$patient_id)) {
    difference <- logcpm[, info$patient_id == patient & info$side == "concave"] -
                  logcpm[, info$patient_id == patient & info$side == "convex"]
    close(difference[expected$gene], expected[[paste0("paired_logCPM_difference.", patient)]],
          paste(lineage, patient), tolerance = 1e-11)
  }
  checked_genes <- checked_genes + nrow(expected)
  cat("PASS pseudobulk refit:", lineage, nrow(expected), "genes\n")
}
verify_model <- function(data_path, output_path, diagnosis = FALSE) {
  data <- read.delim(data_path, check.names = FALSE)
  outputs <- read.delim(output_path, check.names = FALSE)
  data$side <- factor(data$side, levels = c("concave", "convex"))
  if (!diagnosis) data$cohort <- "AIS"
  for (i in seq_len(nrow(outputs))) {
    row <- outputs[i, ]
    subset <- if (diagnosis) data[data$cohort == row$cohort, ] else data
    block <- if (diagnosis) "factor(pair_id)" else "factor(patient_id)"
    formula <- switch(as.character(row$model),
      biology_only = paste("score ~", block, "+ side"),
      biology_plus_library_QC = paste("score ~", block, "+ side + log_library_size"),
      biology_plus_composition_PC1 = paste("score ~", block, "+ side + composition_PC1"),
      QC_only_no_side_term = "score ~ log_library_size + PC1 + PC2", stop("Unknown model"))
    model <- lm(as.formula(formula), subset)
    close(df.residual(model), row$residual_df, paste(data_path, row$model, "df"))
    if (row$model != "QC_only_no_side_term") {
      coefficient <- summary(model)$coefficients["sideconvex", ]
      close(c(-coefficient[1], coefficient[2], coefficient[4]),
        c(row$estimate_concave_minus_convex, row$se, row$p_value), paste(data_path, row$model), tolerance = 1e-10)
    }
  }
  cat("PASS retained score-model refits:", nrow(outputs), "models\n")
}
verify_model("results/tables/bulk/multicellular_axis_model_data.tsv",
             "results/tables/bulk/multicellular_axis_model_sensitivity.tsv")
verify_model("results/tables/disease_control_gse254300_cs/frozen_axis_model_data.tsv",
             "results/tables/disease_control_gse254300_cs/frozen_axis_model_sensitivity.tsv", TRUE)
cat("All five pseudobulk lineages,", checked_genes, "gene fits, donor effects and ten score models passed.\n")
