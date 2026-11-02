#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(edgeR); library(limma) })
source("scripts/lib_common.R")
set.seed(1368770)

raw_root <- raw_data_root()
external_root <- dirname(raw_root)
tx2gene_file <- file.path(external_root, "reference", "gencode_v49",
                          "gencode.v49.tx2gene.tsv")
status <- read.delim("logs/gse311180_requant_status.tsv", stringsAsFactors = FALSE)
stopifnot(nrow(status) == 14L,
          length(unique(status$patient_id)) == 7L,
          all(status$quant_status %in% c("QUANT_PRESENT", "QUANT_PASS")),
          all(file.exists(status$quant_sf)))
tx2gene <- read.delim(tx2gene_file, stringsAsFactors = FALSE)
validate_paired_samples(status, "sample_id", 7L)

read_gene_counts <- function(path) {
  quant <- read.delim(path, stringsAsFactors = FALSE)
  mapped <- merge(quant[, c("Name", "NumReads")], tx2gene,
                  by.x = "Name", by.y = "transcript_id", all.x = TRUE)
  mapping_rate <- mean(!is.na(mapped$gene_name))
  if (mapping_rate < 0.95) stop("Transcript-to-gene mapping below 95% for ", path, ": ", mapping_rate)
  mapped <- mapped[!is.na(mapped$gene_name), ]
  sums <- rowsum(mapped$NumReads, group = mapped$gene_name, reorder = FALSE)
  setNames(as.numeric(sums[, 1]), rownames(sums))
}
vectors <- lapply(status$quant_sf, read_gene_counts)
genes <- sort(unique(unlist(lapply(vectors, names))))
counts <- vapply(vectors, function(x) { out <- setNames(numeric(length(genes)), genes); out[names(x)] <- x; out }, numeric(length(genes)))
rownames(counts) <- genes
colnames(counts) <- status$sample_id
info <- status[match(colnames(counts), status$sample_id), ]
info$side <- factor(info$side, levels = c("convex", "concave"))

keep <- filterByExpr(counts, group = info$side, min.count = 10)
y <- calcNormFactors(DGEList(counts[keep, , drop = FALSE]))
design <- model.matrix(~factor(patient_id) + side, info)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "sideconcave")
de <- topTags(qlf, n = Inf, sort.by = "none")$table
de$gene <- rownames(de)
de$contrast <- "concave_minus_convex"
de$n_paired_patients <- 7L
de$input <- "public_SRA_GENCODEv49_Salmon2.3.4_gene_summed_counts"
write_tsv(de, "results/tables/bulk_gse311180_requant/GSE311180_paired_edgeR_full.tsv.gz")

logcpm <- cpm(y, log = TRUE, prior.count = 2)
write_tsv(data.frame(gene = rownames(logcpm), logcpm),
          "results/tables/bulk_gse311180_requant/GSE311180_TMM_logCPM.tsv.gz")
projection_matrix <- label_blind_logcpm(counts, y$samples$norm.factors)
projection_logcpm <- projection_matrix$logcpm
write_tsv(data.frame(gene = rownames(projection_logcpm), projection_logcpm),
          "results/tables/bulk_gse311180_requant/GSE311180_label_blind_projection_logCPM.tsv.gz")
identity_genes <- intersect(c("CKM", "ACTA1", "MYH1", "MYH2", "DES",
                              "ACAN", "COL2A1", "SOX9"), rownames(logcpm))
identity <- data.frame(sample_id = colnames(logcpm),
                       patient_id = info$patient_id, side = as.character(info$side),
                       t(logcpm[identity_genes, , drop = FALSE]), check.names = FALSE)
muscle_genes <- intersect(c("CKM", "ACTA1", "MYH1", "MYH2", "DES"), names(identity))
cartilage_genes <- intersect(c("ACAN", "COL2A1", "SOX9"), names(identity))
identity$muscle_marker_mean_logCPM <- rowMeans(identity[, muscle_genes, drop = FALSE])
identity$cartilage_marker_mean_logCPM <- rowMeans(identity[, cartilage_genes, drop = FALSE])
identity$muscle_minus_cartilage <- identity$muscle_marker_mean_logCPM -
  identity$cartilage_marker_mean_logCPM
write_tsv(identity, "results/tables/bulk_gse311180_requant/GSE311180_tissue_identity.tsv")
pc <- prcomp(t(logcpm), scale. = FALSE)
read_meta_number <- function(quant_sf, key) {
  meta <- file.path(dirname(quant_sf), "aux_info", "meta_info.json")
  lines <- readLines(meta, warn = FALSE)
  hit <- grep(sprintf('"%s"', key), lines, value = TRUE)
  if (length(hit) != 1L) return(NA_real_)
  as.numeric(sub(".*:[[:space:]]*([0-9.eE+-]+).*", "\\1", hit))
}
qc <- data.frame(
  sample_id = colnames(counts), patient_id = info$patient_id, side = as.character(info$side),
  library_size = colSums(counts), detected_genes = colSums(counts > 0),
  normalization_factor = y$samples$norm.factors,
  read_pairs_processed = vapply(status$quant_sf, read_meta_number, numeric(1), key = "num_processed"),
  percent_mapped = vapply(status$quant_sf, read_meta_number, numeric(1), key = "percent_mapped"),
  PC1 = pc$x[info$sample_id, 1], PC2 = pc$x[info$sample_id, 2]
)
write_tsv(qc, "results/tables/bulk_gse311180_requant/GSE311180_QC.tsv")

signature <- read.delim("data/processed/frozen_multicellular_signature.tsv")
projection <- project_frozen_signature(projection_logcpm, signature)
weights <- projection$used_weights
axis_score <- projection$score
usage <- projection$weights
usage$dataset <- "GSE311180_requant"
usage$detectability_rule <- sprintf("CPM>=%.1f_in_at_least_%d_libraries_label_blind",
                                    projection_matrix$min_cpm,
                                    projection_matrix$min_samples)
write_tsv(usage, "results/tables/bulk_gse311180_requant/frozen_axis_projection_genes.tsv")
axis <- data.frame(
  sample_id = names(axis_score), score = as.numeric(axis_score),
  patient_id = info$patient_id[match(names(axis_score), info$sample_id)],
  side = as.character(info$side[match(names(axis_score), info$sample_id)])
)
axis_eff <- merge(axis[axis$side == "concave", ], axis[axis$side == "convex", ],
                  by = "patient_id", suffixes = c(".concave", ".convex"))
axis_eff$paired_difference <- axis_eff$score.concave - axis_eff$score.convex
grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(axis_eff))))
perm <- as.numeric(grid %*% axis_eff$paired_difference / nrow(axis_eff))
axis_test <- data.frame(
  n_pairs = nrow(axis_eff), mean_paired_difference = mean(axis_eff$paired_difference),
  se = sd(axis_eff$paired_difference) / sqrt(nrow(axis_eff)),
  t_test_p = t.test(axis_eff$paired_difference)$p.value,
  exact_sign_flip_p = mean(abs(perm) >= abs(mean(axis_eff$paired_difference)) - 1e-12),
  positive_pairs = sum(axis_eff$paired_difference > 0)
)
write_tsv(axis, "results/tables/bulk_gse311180_requant/frozen_axis_scores.tsv")
write_tsv(axis_eff, "results/tables/bulk_gse311180_requant/frozen_axis_paired_effects.tsv")
write_tsv(axis_test, "results/tables/bulk_gse311180_requant/frozen_axis_test.tsv")

# Null signatures preserve signature size, expression decile and the empirical
# weight distribution. They test whether the selected cross-lineage genes are
# more coherent than equivalently abundant random genes in this new cohort.
mean_expr <- rowMeans(projection_logcpm)
expr_decile <- cut(rank(mean_expr, ties.method = "average"),
                   breaks = quantile(rank(mean_expr, ties.method = "average"),
                                     probs = seq(0, 1, .1)),
                   include.lowest = TRUE, labels = FALSE)
names(expr_decile) <- names(mean_expr)
weight_decile <- expr_decile[weights$gene]
random_axis_means <- vapply(seq_len(1000), function(iteration) {
  random_genes <- sample_matched_genes_without_replacement(
    weights$gene, weight_decile, expr_decile, weights$gene)
  random_z <- t(scale(t(projection_logcpm[random_genes, , drop = FALSE])))
  random_z[!is.finite(random_z)] <- 0
  random_weights <- sample(weights$weight)
  random_score <- colSums(random_z * random_weights) / sum(abs(random_weights))
  random_df <- data.frame(sample_id = names(random_score), score = random_score,
                          patient_id = info$patient_id[match(names(random_score), info$sample_id)],
                          side = as.character(info$side[match(names(random_score), info$sample_id)]))
  random_eff <- merge(random_df[random_df$side == "concave", ],
                      random_df[random_df$side == "convex", ], by = "patient_id",
                      suffixes = c(".concave", ".convex"))
  mean(random_eff$score.concave - random_eff$score.convex)
}, numeric(1))
random_summary <- data.frame(
  observed_mean = mean(axis_eff$paired_difference), iterations = length(random_axis_means),
  expression_decile_and_size_matched_empirical_p =
    (1 + sum(abs(random_axis_means) >= abs(mean(axis_eff$paired_difference)))) /
    (1 + length(random_axis_means))
)
write_tsv(data.frame(iteration = seq_along(random_axis_means),
                     random_mean_paired_difference = random_axis_means,
                     n_unique_genes = nrow(weights),
                     target_unique_genes = nrow(weights)),
          "results/tables/bulk_gse311180_requant/frozen_axis_random_signature_null.tsv.gz")
write_tsv(random_summary,
          "results/tables/bulk_gse311180_requant/frozen_axis_random_signature_test.tsv")

# Repeat the frozen-axis analysis after removing genes in the top one per cent
# of cohort mean expression.
high_cut <- quantile(mean_expr, .99)
trimmed_weights <- weights[mean_expr[weights$gene] < high_cut, ]
trimmed_z <- t(scale(t(projection_logcpm[trimmed_weights$gene, , drop = FALSE])))
trimmed_z[!is.finite(trimmed_z)] <- 0
trimmed_score <- colSums(trimmed_z * trimmed_weights$weight) /
  sum(abs(trimmed_weights$weight))
trimmed <- data.frame(sample_id = names(trimmed_score), score = trimmed_score,
                      patient_id = info$patient_id[match(names(trimmed_score), info$sample_id)],
                      side = as.character(info$side[match(names(trimmed_score), info$sample_id)]))
trimmed_eff <- merge(trimmed[trimmed$side == "concave", ],
                     trimmed[trimmed$side == "convex", ], by = "patient_id",
                     suffixes = c(".concave", ".convex"))
trimmed_eff$paired_difference <- trimmed_eff$score.concave - trimmed_eff$score.convex
write_tsv(data.frame(n_signature_genes = nrow(weights),
                     n_retained_genes = nrow(trimmed_weights),
                     mean_paired_difference = mean(trimmed_eff$paired_difference),
                     positive_pairs = sum(trimmed_eff$paired_difference > 0),
                     t_test_p = t.test(trimmed_eff$paired_difference)$p.value),
          "results/tables/bulk_gse311180_requant/frozen_axis_remove_top1pct.tsv")

frozen <- readRDS("data/processed/frozen_pathway_definitions.rds")
path_scores <- gene_set_zscore(logcpm, frozen$pathways, min_genes = 5)
path_long <- do.call(rbind, lapply(rownames(path_scores), function(gs) data.frame(
  gene_set = gs, sample_id = colnames(path_scores), score = path_scores[gs, ],
  patient_id = info$patient_id[match(colnames(path_scores), info$sample_id)],
  side = as.character(info$side[match(colnames(path_scores), info$sample_id)])
)))
path_eff <- merge(path_long[path_long$side == "concave", ], path_long[path_long$side == "convex", ],
                  by = c("gene_set", "patient_id"), suffixes = c(".concave", ".convex"))
path_eff$paired_difference <- path_eff$score.concave - path_eff$score.convex
path_tests <- do.call(rbind, lapply(split(path_eff, path_eff$gene_set), function(z) {
  obs <- mean(z$paired_difference)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(z))))
  null <- as.numeric(signs %*% z$paired_difference / nrow(z))
  data.frame(gene_set = z$gene_set[1], n_pairs = nrow(z), mean_paired_difference = obs,
             positive_pairs = sum(z$paired_difference > 0),
             t_test_p = t.test(z$paired_difference)$p.value,
             exact_sign_flip_p = mean(abs(null) >= abs(obs) - 1e-12))
}))
path_tests$FDR_t_test <- p.adjust(path_tests$t_test_p, "BH")
path_tests$FDR_exact <- p.adjust(path_tests$exact_sign_flip_p, "BH")
write_tsv(path_long, "results/tables/bulk_gse311180_requant/frozen_pathway_scores.tsv")
write_tsv(path_eff, "results/tables/bulk_gse311180_requant/frozen_pathway_paired_effects.tsv")
write_tsv(path_tests, "results/tables/bulk_gse311180_requant/frozen_pathway_paired_tests.tsv")

loo <- do.call(rbind, lapply(unique(axis_eff$patient_id), function(left_out) {
  z <- axis_eff[axis_eff$patient_id != left_out, ]
  data.frame(left_out_patient = left_out, n_pairs = nrow(z),
             mean_paired_difference = mean(z$paired_difference),
             positive_pairs = sum(z$paired_difference > 0),
             t_test_p = t.test(z$paired_difference)$p.value)
}))
write_tsv(loo, "results/tables/bulk_gse311180_requant/frozen_axis_leave_one_patient_out.tsv")

# Cross-cohort replication and patient-level combination: 5 GSE254300 pairs + 7 GSE311180 pairs.
old <- read.delim(gzfile(analysis_input("results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz")))
replication <- merge(old[, c("gene", "logFC", "PValue", "FDR")],
                     de[, c("gene", "logFC", "PValue", "FDR")], by = "gene",
                     suffixes = c(".GSE254300", ".GSE311180_requant"))
replication$sign_concordant <- sign(replication$logFC.GSE254300) == sign(replication$logFC.GSE311180_requant)
z1 <- sign(replication$logFC.GSE254300) * qnorm(pmax(replication$PValue.GSE254300 / 2, 1e-300), lower.tail = FALSE)
z2 <- sign(replication$logFC.GSE311180_requant) * qnorm(pmax(replication$PValue.GSE311180_requant / 2, 1e-300), lower.tail = FALSE)
replication$meta_z <- (sqrt(5) * z1 + sqrt(7) * z2) / sqrt(12)
replication$meta_p <- 2 * pnorm(abs(replication$meta_z), lower.tail = FALSE)
replication$meta_FDR <- p.adjust(replication$meta_p, "BH")
write_tsv(replication, "results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz")
write_tsv(data.frame(
  metric = c("spearman_logFC_all", "sign_concordance_all", "sign_concordance_GSE254300_FDR05",
             "sign_concordance_GSE311180_FDR05"),
  value = c(cor(replication$logFC.GSE254300, replication$logFC.GSE311180_requant, method = "spearman"),
            mean(replication$sign_concordant),
            mean(replication$sign_concordant[replication$FDR.GSE254300 < .05], na.rm = TRUE),
            mean(replication$sign_concordant[replication$FDR.GSE311180_requant < .05], na.rm = TRUE)),
  n = c(nrow(replication), nrow(replication), sum(replication$FDR.GSE254300 < .05),
        sum(replication$FDR.GSE311180_requant < .05))
), "results/tables/meta_validation_requant/gene_replication_summary.tsv")

old_axis <- read.delim("results/tables/bulk/multicellular_axis_paired_effects.tsv")
combined_axis <- rbind(
  data.frame(cohort = "GSE254300", patient_id = old_axis$patient_id,
             paired_difference = old_axis$paired_difference),
  data.frame(cohort = "GSE311180_requant", patient_id = axis_eff$patient_id,
             paired_difference = axis_eff$paired_difference)
)
# Each constituent score is already constructed after gene-wise z-scoring
# within its own cohort. The paired differences are therefore on the same
# standardized score scale; no effect-level re-standardization is performed.
combined_axis$score_scale <- "cohort_internal_gene_z_scores"
grid12 <- as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(combined_axis))))
null12 <- as.numeric(grid12 %*% combined_axis$paired_difference / nrow(combined_axis))
combined_axis_test <- data.frame(
  n_patients = nrow(combined_axis), mean_paired_difference = mean(combined_axis$paired_difference),
  se = sd(combined_axis$paired_difference) / sqrt(nrow(combined_axis)),
  t_test_p = t.test(combined_axis$paired_difference)$p.value,
  exact_sign_flip_p = mean(abs(null12) >= abs(mean(combined_axis$paired_difference)) - 1e-12),
  positive_pairs = sum(combined_axis$paired_difference > 0)
)
write_tsv(combined_axis, "results/tables/meta_validation_requant/frozen_axis_12_patient_effects.tsv")
write_tsv(combined_axis_test, "results/tables/meta_validation_requant/frozen_axis_12_patient_test.tsv")

capture_versions(c("edgeR", "limma"), "logs/gse311180_requant_R_versions.tsv")
message("GSE311180 raw requantification validation complete for seven patient pairs.")
