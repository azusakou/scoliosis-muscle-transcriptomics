#!/usr/bin/env Rscript
# Recompute downstream figure inputs from published retained records. This does
# not rerun cell calling, SoupX, annotation, embeddings, LIANA or RNA alignment.
# Reference tables remain immutable; all derived checks are written under logs.
source("scripts/lib_common.R")
out <- "logs/retained_figure_inputs"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
rd <- function(path) read.delim(if (grepl("\\.gz$", path)) gzfile(path) else path,
                              check.names = FALSE, stringsAsFactors = FALSE)
sha <- function(path) {
  command <- if (nzchar(Sys.which("sha256sum"))) "sha256sum" else "shasum"
  options <- if (command == "shasum") c("-a", "256") else character()
  value <- system2(command, c(options, shQuote(path)), stdout = TRUE)
  if (!is.null(attr(value, "status"))) stop("Cannot hash figure input: ", path)
  substr(value[[1]], 1, 64)
}
checks <- list()
key <- function(x, columns) do.call(paste, c(x[columns], sep = "\034"))
same <- function(actual, expected, keys, columns, name, exact = character()) {
  a <- key(actual, keys); b <- key(expected, keys)
  if (anyDuplicated(a) || anyDuplicated(b) || !setequal(a, b))
    stop(name, ": changed, missing or duplicate source keys")
  actual <- actual[match(b, a), , drop = FALSE]
  for (column in columns) {
    x <- actual[[column]]; y <- expected[[column]]
    if (!identical(is.na(x), is.na(y))) stop(name, ": changed missingness in ", column)
    keep <- !is.na(x)
    if (is.numeric(x) && is.numeric(y) && !column %in% exact) {
      if (any(!is.finite(x[keep])) || any(!is.finite(y[keep])) ||
          any(abs(x[keep] - y[keep]) > 1e-12 + 1e-10 * abs(y[keep])))
        stop(name, ": numerical mismatch in ", column)
    } else if (!identical(as.character(x[keep]), as.character(y[keep]))) {
      stop(name, ": exact count/label mismatch in ", column)
    }
  }
  checks[[length(checks) + 1L]] <<- data.frame(check = name, rows = nrow(expected),
    columns_checked = paste(columns, collapse = ","), status = "PASS",
    scope = "Recomputed from retained inputs; upstream algorithms not rerun")
  invisible(actual)
}
same_multiset <- function(actual, expected, keys, columns, name) {
  # Some legacy LIANA complex records share gene-level keys. Preserve every
  # original merged row, sorting within keys to compare the full multiset.
  repeated <- c(actual = sum(duplicated(key(actual, keys))),
                expected = sum(duplicated(key(expected, keys))))
  if (nrow(actual) != nrow(expected) || repeated[1] != repeated[2])
    stop(name, ": changed legacy row multiplicity")
  annotate <- function(x) {
    x <- x[do.call(order, unname(as.list(x[c(keys, columns)]))), , drop = FALSE]
    x$.occurrence <- ave(seq_len(nrow(x)), key(x, keys), FUN = seq_along)
    x
  }
  same(annotate(actual), annotate(expected), c(keys, ".occurrence"), columns, name)
  write_tsv(data.frame(check = name, rows = nrow(actual), repeated_gene_key_rows = repeated[1],
    scope = "Legacy complex-to-gene merge multiplicity preserved; rows are not independent evidence"),
    file.path(out, "communication_legacy_row_multiplicity.tsv"))
}
short <- c("fibro-adipogenic progenitor" = "FAP", "endothelial cell" = "Endothelial",
  "macrophage/monocyte" = "Mac/Mono", "muscle stem/progenitor cell" = "MuSC",
  "smooth-muscle-like cell" = "SMC-like")

# Cell metadata are sufficient for composition, with all verified zeros retained.
metadata <- rd("results/tables/annotation/annotated_cell_metadata.tsv.gz")
stopifnot(nrow(metadata) == 14846L, !anyDuplicated(metadata$cell_id),
          all(metadata$passes_qc), all(metadata$scDblFinder.class == "singlet"))
counts <- composition_from_cells(metadata)
same(counts, rd("results/tables/composition/cell_type_counts_and_proportions.tsv"),
  c("sample_id", "cell_type"), c("patient_id", "side", "Freq", "sample_total", "proportion"),
  "composition_from_retained_cell_labels", c("Freq", "sample_total"))
paired <- reshape(counts[, c("patient_id", "side", "cell_type", "proportion")],
  idvar = c("patient_id", "cell_type"), timevar = "side", direction = "wide")
paired$concave_minus_convex <- paired$proportion.concave - paired$proportion.convex
same(paired, rd("results/tables/composition/paired_composition_effects.tsv"),
  c("patient_id", "cell_type"), c("proportion.concave", "proportion.convex", "concave_minus_convex"),
  "paired_composition_including_zero_counts")
write_tsv(counts, file.path(out, "composition_from_cells.tsv"))
write_tsv(paired, file.path(out, "paired_composition.tsv"))

# The two published coordinate tables must describe the same retained cells.
umap <- rd("results/tables/annotation/annotated_umap.tsv")
qc_umap <- rd("results/tables/qc/scrna_unintegrated_umap.tsv")
same(umap, metadata, "cell_id", c("patient_id", "side", "cell_type"), "atlas_cell_label_join")
same(qc_umap, umap, "cell_id", c("patient_id", "side", "umapunintegrated_1", "umapunintegrated_2"),
  "same_archived_coordinates_in_QC_and_atlas")
qc <- rd("results/tables/qc/scrna_qc_summary.tsv")
totals <- aggregate(Freq ~ sample_id + patient_id + side, counts, sum)
totals$final_singlets <- totals$Freq
same(totals, qc, "sample_id", c("patient_id", "side", "final_singlets"),
  "retained_library_counts", "final_singlets")

# Re-pair the supplied gene-set scores. No GMT or new gene-set scoring is implied.
scores <- rd("results/tables/programs/scrna_pathway_scores.tsv")
keys <- c("cell_type", "gene_set", "patient_id")
effects <- merge(scores[scores$side == "concave", ], scores[scores$side == "convex", ],
  by = keys, suffixes = c(".concave", ".convex"))
stopifnot(!anyDuplicated(key(effects, keys)))
effects$paired_difference <- effects$score.concave - effects$score.convex
same(effects, rd("results/tables/programs/scrna_pathway_paired_effects_by_donor.tsv"),
  keys, c("score.concave", "score.convex", "paired_difference"), "programme_paired_scores")
summary <- aggregate(paired_difference ~ cell_type + gene_set, effects,
  function(x) c(mean = mean(x), concordant = as.numeric(all(x > 0) || all(x < 0))))
summary <- cbind(summary[1:2], as.data.frame(summary$paired_difference))
same(summary, rd("results/tables/programs/scrna_pathway_effect_summary.tsv"),
  c("cell_type", "gene_set"), c("mean", "concordant"), "programme_mean_and_direction", "concordant")
core <- c("HALLMARK_MYOGENESIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_ANGIOGENESIS",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
display <- effects[effects$gene_set %in% core, ]
stopifnot(nrow(display) == 120L)
same(display, rd("results/figures/source_data/Figure_3_donor_pathways.tsv"), keys,
  "paired_difference", "displayed12_programme_donor_effects")
concordance <- aggregate(concordant ~ cell_type, summary[summary$gene_set %in% core, ], mean)
concordance$concordant_pct <- 100 * concordance$concordant
same(concordance, rd("results/figures/source_data/Figure_3_concordance.tsv"), "cell_type",
  c("concordant", "concordant_pct"), "displayed12_programme_concordance")
write_tsv(effects, file.path(out, "programme_paired_scores.tsv"))
write_tsv(summary, file.path(out, "programme_direction_summary.tsv"))
write_tsv(concordance, file.path(out, "programme_core_concordance.tsv"))

# Only means of the retained component-level R2 metrics are reconstructed.
embedding <- rd("results/tables/robustness/integration_embedding_sensitivity.tsv")
embedding_summary <- aggregate(cbind(side_semipartial_r2, sample_r2) ~ method, embedding, mean)
same(embedding_summary, rd("results/tables/robustness/integration_embedding_sensitivity_summary.tsv"),
  "method", c("side_semipartial_r2", "sample_r2"), "embedding_component_summary_only")
write_tsv(embedding_summary, file.path(out, "embedding_component_means.tsv"))

# Reapply the exact original script50 postprocessing to complete retained LIANA
# method outputs, without invoking LIANA or manufacturing cell expression.
all <- rd("results/tables/communication/LIANA_NATMI_CellPhoneDB_all.tsv.gz")
stopifnot(nrow(all) == 144872L, all(all$method %in% c("natmi", "cellphonedb")))
expected_score <- ifelse(all$method == "natmi", all$edge_specificity, all$lr.mean)
expected_support <- ifelse(all$method == "natmi", TRUE, all$pvalue <= .05)
same(data.frame(row = seq_len(nrow(all)), method_score = expected_score),
  data.frame(row = seq_len(nrow(all)), method_score = all$method_score), "row",
  "method_score", "LIANA_native_score_columns")
# The original native caches were independently matched against all 144,872
# published rows/flags. Their 37 actual boundary P values are supplied at 17
# significant digits, restoring the original <=.05 decision without rounding.
boundary <- all$method == "cellphonedb" & !is.na(all$pvalue) & all$pvalue == .05
native_keys <- c("method", "sample_id", "source", "target", "ligand.complex", "ligand", "receptor.complex", "receptor")
if ("native_pvalue_17g" %in% names(all)) {
  # A future raw script50 output carries its complete native P precision. It
  # must not reuse the old published-table precision map after source changes.
  cp <- all$method == "cellphonedb"
  native_p <- as.numeric(all$native_pvalue_17g[cp])
  stopifnot(all(is.finite(native_p)), all(native_p >= 0 & native_p <= 1),
            all(abs(native_p - all$pvalue[cp]) <= 1e-12))
  expected_support[cp] <- native_p <= .05
  precision <- all[boundary, c(native_keys, "patient_id", "side", "pvalue", "native_pvalue_17g", "method_supported")]
} else {
precision <- read.delim("analysis/communication_native_precision_boundaries.tsv", check.names = FALSE,
  colClasses = "character", stringsAsFactors = FALSE)
native_sources <- rd("analysis/communication_native_precision_sources.tsv")
stopifnot(sum(boundary) == 37L, nrow(precision) == 37L,
          identical(as.integer(precision$published_row), which(boundary)),
          identical(key(precision, native_keys), key(all[boundary, ], native_keys)),
          !anyDuplicated(native_sources$sample_id), nrow(native_sources) == 4L,
          all(native_sources$published_rows_matched == nrow(all)),
          all(precision$published_all_table_sha256 == sha("results/tables/communication/LIANA_NATMI_CellPhoneDB_all.tsv.gz")),
          all(native_sources$published_all_table_sha256 == precision$published_all_table_sha256[[1]]),
          all(native_sources$exporter_sha256 == sha("scripts/export_communication_precision.R")))
source_index <- match(precision$sample_id, native_sources$sample_id)
stopifnot(!anyNA(source_index),
          identical(precision$cache_basename, native_sources$cache_basename[source_index]),
          identical(precision$cache_sha256, native_sources$cache_sha256[source_index]))
native_p <- as.numeric(precision$native_pvalue_17g)
stopifnot(all(is.finite(native_p)), all(native_p > .05),
          all(native_p <= .05 + 1e-14), all(precision$method_supported == "FALSE"))
expected_support[boundary] <- native_p <= .05
}
same(data.frame(row = seq_len(nrow(all)), method_supported = expected_support),
  data.frame(row = seq_len(nrow(all)), method_supported = all$method_supported), "row",
  "method_supported", "LIANA_native_precision_support_for_all_rows")
write_tsv(precision, file.path(out, "communication_serialized_P_boundary.tsv"))
z <- all[all$ligand.prop >= .1 & all$receptor.prop >= .1 & all$source != all$target & all$method_supported, ]
keys <- c("method", "source", "target", "ligand", "receptor", "patient_id")
paired_edges <- merge(z[z$side == "concave", c(keys, "method_score")],
  z[z$side == "convex", c(keys, "method_score")], by = keys, suffixes = c(".concave", ".convex"))
paired_edges$paired_score_difference <- paired_edges$method_score.concave - paired_edges$method_score.convex
same_multiset(paired_edges, rd("results/tables/communication/paired_method_edges_by_donor.tsv.gz"), keys,
  c("method_score.concave", "method_score.convex", "paired_score_difference"), "communication_per_method_donor_pairs")
method <- aggregate(paired_score_difference ~ method + source + target + ligand + receptor, paired_edges,
  function(x) c(mean = mean(x), donor_concordant = as.numeric(length(x) == 2 && (all(x > 0) || all(x < 0))), n = length(x)))
method <- cbind(method[1:5], as.data.frame(method$paired_score_difference))
cross <- merge(method[method$method == "natmi", ], method[method$method == "cellphonedb", ],
  by = c("source", "target", "ligand", "receptor"), suffixes = c(".natmi", ".cellphonedb"))
cross$method_direction_concordant <- sign(cross$mean.natmi) == sign(cross$mean.cellphonedb)
pb <- rd("results/tables/pseudobulk/scrna_donor_specific_effects.tsv.gz")
pb_keys <- key(pb, c("gene", "cell_type", "patient_id"))
stopifnot(!anyDuplicated(pb_keys))
for (pid in unique(pb$patient_id)) {
  ligand <- paste(cross$ligand, cross$source, pid, sep = "\034")
  receptor <- paste(cross$receptor, cross$target, pid, sep = "\034")
  cross[[paste0("pseudobulk_joint_effect_", pid)]] <-
    pb$paired_logCPM_difference[match(ligand, pb_keys)] +
    pb$paired_logCPM_difference[match(receptor, pb_keys)]
}
joint <- grep("^pseudobulk_joint_effect", names(cross), value = TRUE)
cross$pseudobulk_donor_concordant <- apply(cross[, joint, drop = FALSE], 1,
  function(x) all(is.finite(x)) && (all(x > 0) || all(x < 0)))
cross$pseudobulk_matches_methods <- cross$pseudobulk_donor_concordant &
  sign(rowMeans(cross[, joint, drop = FALSE], na.rm = TRUE)) == sign(cross$mean.natmi)
cross$retained <- cross$donor_concordant.natmi == 1 & cross$donor_concordant.cellphonedb == 1 &
  cross$method_direction_concordant & cross$pseudobulk_matches_methods
cross <- cross[order(cross$retained, abs(cross$mean.natmi) + abs(cross$mean.cellphonedb), decreasing = TRUE), ]
columns <- c("mean.natmi", "donor_concordant.natmi", "n.natmi", "mean.cellphonedb", "donor_concordant.cellphonedb",
  "n.cellphonedb", "method_direction_concordant", joint, "pseudobulk_donor_concordant", "pseudobulk_matches_methods", "retained")
same(cross, rd("results/tables/communication/cross_method_donor_aware_edges.tsv"),
  c("source", "target", "ligand", "receptor"), columns, "communication_original_filter_reconstruction",
  c("donor_concordant.natmi", "n.natmi", "donor_concordant.cellphonedb", "n.cellphonedb"))
stopifnot(nrow(cross) == 7516L, sum(cross$retained) == 286L)
selected <- cross[cross$retained, ]; selected$combined <- sign(selected$mean.natmi) *
  (abs(selected$mean.natmi) + abs(selected$mean.cellphonedb))
selected <- head(selected[order(abs(selected$combined), decreasing = TRUE), ], 10)
same(selected, rd("results/figures/source_data/Figure_5_communication.tsv"),
  c("source", "target", "ligand", "receptor"), c("mean.natmi", "mean.cellphonedb", "combined"),
  "displayed10_communication_candidate_selection")
write_tsv(cross, file.path(out, "communication_cross_method_reconstruction.tsv"))
write_tsv(selected, file.path(out, "communication_selected10.tsv"))

# Target-set panel selection follows the last submitted source producer84:
# two largest absolute concordant means per lineage, from supplied target scores.
tf <- rd("results/tables/programs/scrna_gtrd_tf_target_scores_top_variable.tsv")
tf_effect <- merge(tf[tf$side == "concave", ], tf[tf$side == "convex", ],
  by = c("cell_type", "tf_target_set", "patient_id"), suffixes = c(".concave", ".convex"))
tf_effect$effect <- tf_effect$score.concave - tf_effect$score.convex
tf_summary <- aggregate(effect ~ cell_type + tf_target_set, tf_effect,
  function(x) c(mean = mean(x), min = min(x), max = max(x), concordant = as.numeric(all(x > 0) || all(x < 0))))
tf_summary <- cbind(tf_summary[1:2], as.data.frame(tf_summary$effect))
tf_summary <- tf_summary[tf_summary$concordant == 1, ]
tf_selected <- do.call(rbind, lapply(split(tf_summary, tf_summary$cell_type),
  function(x) head(x[order(abs(x$mean), decreasing = TRUE), ], 2)))
same(tf_selected, rd("results/figures/source_data/Figure_5_TF.tsv"), c("cell_type", "tf_target_set"),
  c("mean", "min", "max", "concordant"), "displayed_target_set_selection", "concordant")
write_tsv(tf_selected, file.path(out, "target_set_selected.tsv"))
write_tsv(do.call(rbind, checks), file.path(out, "checks.tsv"))
writeLines(c("Retained downstream figure-input reconstruction PASS.",
  "Original reference tables were read and compared, never overwritten.",
  "Method-support decisions use round-trip native P precision, preserving the original <=.05 rule for every method row.",
  "No new EmptyDrops, SoupX, annotation, gene-set scoring, PCA/UMAP/RPCA/Harmony, LIANA, Salmon or STAR run."),
  file.path(out, "SCOPE.txt"))
cat("Retained figure-input checks PASS:", length(checks), "comparisons. See", out, "\n")
