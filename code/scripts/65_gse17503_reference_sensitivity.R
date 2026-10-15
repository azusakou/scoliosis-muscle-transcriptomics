#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))
source("scripts/lib_common.R")

raw_root <- raw_data_root()
raw_dir <- file.path(raw_root, "reference_gse17503")
matrix_file <- file.path(raw_dir, "GSE17503_series_matrix.txt.gz")
annotation_file <- file.path(raw_dir, "GSE17503_non-normalized_data.txt.gz")
stopifnot(file.exists(matrix_file), file.exists(annotation_file))

matrix_lines <- readLines(gzfile(matrix_file), warn = FALSE)
begin <- grep("^!series_matrix_table_begin", matrix_lines)
end <- grep("^!series_matrix_table_end", matrix_lines)
stopifnot(length(begin) == 1L, length(end) == 1L, end > begin)
con <- textConnection(matrix_lines[(begin + 1):(end - 1)])
on.exit(close(con), add = TRUE)
probe_table <- read.delim(con, check.names = FALSE, quote = "\"",
                          stringsAsFactors = FALSE)
rownames(probe_table) <- probe_table$ID_REF
probe_matrix <- as.matrix(probe_table[, -1, drop = FALSE])
storage.mode(probe_matrix) <- "double"

get_meta <- function(key) {
  line <- matrix_lines[grep(paste0("^", key, "\\t"), matrix_lines)[1]]
  values <- strsplit(line, "\\t", fixed = FALSE)[[1]][-1]
  gsub('^"|"$', "", values)
}
sample_titles <- get_meta("!Sample_title")
sample_accessions <- get_meta("!Sample_geo_accession")
stopifnot(identical(sample_accessions, colnames(probe_matrix)), length(sample_titles) == 10L)
sample_info <- data.frame(
  sample_id = sample_accessions,
  title = sample_titles,
  patient_id = sub(" .*", "", sample_titles),
  material = ifelse(grepl("vitro", sample_titles, ignore.case = TRUE), "culture", "biopsy"),
  stringsAsFactors = FALSE
)
stopifnot(length(unique(sample_info$patient_id)) == 5L,
          all(table(sample_info$patient_id) == 2L),
          all(table(sample_info$material) == 5L))

annotation <- read.delim(gzfile(annotation_file), check.names = FALSE,
                         stringsAsFactors = FALSE)
stopifnot(all(c("ID_REF", "TRANSCRIPT") %in% names(annotation)))
# The legacy Illumina export is shifted relative to its final annotation
# headers: the gene symbol is stored in TRANSCRIPT (verified against TargetID),
# whereas the field labelled SYMBOL contains a RefSeq accession.
mapping <- unique(annotation[, c("ID_REF", "TRANSCRIPT")])
names(mapping)[2] <- "gene_symbol"
mapping <- mapping[nzchar(mapping$gene_symbol) & mapping$ID_REF %in% rownames(probe_matrix), ]

# The deposited series matrix contains quantile-normalised signal intensities,
# not log values. Negative background-corrected values are truncated at zero
# before log2(x + 1). For duplicated symbols, retain the probe with the largest
# across-sample mean, a prespecified and reproducible specificity rule.
log_expr <- log2(pmax(probe_matrix[mapping$ID_REF, , drop = FALSE], 0) + 1)
probe_mean <- rowMeans(log_expr)
mapping$probe_mean <- probe_mean[mapping$ID_REF]
mapping <- mapping[order(mapping$gene_symbol, -mapping$probe_mean, mapping$ID_REF), ]
mapping <- mapping[!duplicated(mapping$gene_symbol), ]
gene_expr <- log_expr[mapping$ID_REF, , drop = FALSE]
rownames(gene_expr) <- mapping$gene_symbol

sample_info$material <- factor(sample_info$material, levels = c("culture", "biopsy"))
design <- model.matrix(~factor(patient_id) + material, sample_info)
fit <- eBayes(lmFit(gene_expr, design), robust = TRUE)
de <- topTable(fit, coef = "materialbiopsy", number = Inf, sort.by = "none")
de$gene <- rownames(de)
de$contrast <- "biopsy_minus_patient_matched_culture"
de$n_paired_patients <- 5L
write_tsv(de, "results/tables/reference_gse17503/GSE17503_paired_limma_full.tsv.gz")
write_tsv(sample_info, "metadata/audit/gse17503_samples.tsv")

marker_reference <- read.delim("results/tables/bulk/deconvolution_marker_reference.tsv")
marker_sets <- split(marker_reference$gene, marker_reference$cell_type)
z <- t(scale(t(gene_expr)))
z[!is.finite(z)] <- 0
scores <- do.call(rbind, lapply(names(marker_sets), function(cell_type) {
  genes <- intersect(unique(marker_sets[[cell_type]]), rownames(z))
  if (length(genes) < 5) return(NULL)
  data.frame(cell_type = cell_type, n_genes = length(genes),
             sample_id = colnames(z), score = colMeans(z[genes, , drop = FALSE]))
}))
stopifnot(is.data.frame(scores), "sample_id" %in% names(scores))
scores$patient_id <- sample_info$patient_id[match(scores$sample_id, sample_info$sample_id)]
scores$title <- sample_info$title[match(scores$sample_id, sample_info$sample_id)]
scores$material <- as.character(sample_info$material[match(scores$sample_id, sample_info$sample_id)])
effects <- merge(scores[scores$material == "biopsy", ],
                 scores[scores$material == "culture", ],
                 by = c("cell_type", "n_genes", "patient_id"),
                 suffixes = c(".biopsy", ".culture"))
effects$paired_difference <- effects$score.biopsy - effects$score.culture
tests <- do.call(rbind, lapply(split(effects, effects$cell_type), function(x) {
  observed <- mean(x$paired_difference)
  grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), nrow(x))))
  null <- as.numeric(grid %*% x$paired_difference / nrow(x))
  data.frame(cell_type = x$cell_type[1], n_genes = x$n_genes[1], n_pairs = nrow(x),
             mean_biopsy_minus_culture = observed,
             positive_pairs = sum(x$paired_difference > 0),
             t_test_p = t.test(x$paired_difference)$p.value,
             exact_sign_flip_p = mean(abs(null) >= abs(observed) - 1e-12))
}))
tests$FDR_t_test <- p.adjust(tests$t_test_p, "BH")
tests$FDR_exact <- p.adjust(tests$exact_sign_flip_p, "BH")
write_tsv(scores, "results/tables/reference_gse17503/cell_type_marker_scores.tsv")
write_tsv(effects, "results/tables/reference_gse17503/cell_type_marker_paired_effects.tsv")
write_tsv(tests, "results/tables/reference_gse17503/cell_type_marker_paired_tests.tsv")

# This projection is a culture-bias diagnostic only. It must not be interpreted
# as concave-versus-convex replication because GSE17503 has no side labels.
signature <- read.delim("data/processed/frozen_multicellular_signature.tsv")
weights <- aggregate(weight ~ gene, signature, mean)
weights <- weights[weights$gene %in% rownames(z), ]
signature_score <- colSums(z[weights$gene, , drop = FALSE] * weights$weight) /
  sum(abs(weights$weight))
signature_scores <- merge(data.frame(sample_id = names(signature_score),
                                     score = as.numeric(signature_score)),
                          sample_info, by = "sample_id")
signature_effects <- merge(signature_scores[signature_scores$material == "biopsy", ],
                           signature_scores[signature_scores$material == "culture", ],
                           by = "patient_id", suffixes = c(".biopsy", ".culture"))
signature_effects$paired_difference <- signature_effects$score.biopsy -
  signature_effects$score.culture
write_tsv(signature_scores, "results/tables/reference_gse17503/frozen_axis_scores_context_only.tsv")
write_tsv(signature_effects,
          "results/tables/reference_gse17503/frozen_axis_paired_effects_context_only.tsv")

capture_versions(c("limma"), "logs/gse17503_reference_R_versions.tsv")
message("GSE17503 five-patient biopsy/culture reference sensitivity complete.")
