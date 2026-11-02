#!/usr/bin/env Rscript
source("scripts/lib_common.R")

# This step reads discovery outputs only and does not use bulk expression or
# side labels in its calculations. The current script order does not establish
# historical validation-blinded selection or a prospective freeze.
# Manifest input-use flags describe this function, not historical access order.
de <- read.delim(gzfile("results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz"),
                 check.names = FALSE)
eligible <- de$donor_sign_concordant %in% TRUE &
  !grepl("^(MT-|RPS|RPL|HBA|HBB)", de$gene)
de <- de[eligible, , drop = FALSE]

signature <- do.call(rbind, lapply(split(de, de$cell_type), function(z) {
  z <- z[order(abs(z$mean_donor_effect), decreasing = TRUE), , drop = FALSE]
  z <- head(z, 30L)
  z$rank_within_lineage <- seq_len(nrow(z))
  z
}))
signature$weight <- sign(signature$mean_donor_effect)
signature$selection_source <- "PRJNA722100_discovery_only"
signature$selection_rule <- "top30_absolute_mean_donor_effect_same_direction_in_both_donors"
signature <- signature[, c("gene", "cell_type", "rank_within_lineage",
                           "mean_donor_effect", "weight", "donor_sign_concordant",
                           "selection_source", "selection_rule")]

stopifnot(nrow(signature) == 150L,
          length(unique(signature$cell_type)) == 5L,
          all(table(signature$cell_type) == 30L),
          all(signature$weight %in% c(-1, 1)))

weights <- aggregate(weight ~ gene, signature, mean)
weights$n_lineages <- as.integer(table(signature$gene)[weights$gene])
weights$selection_source <- "PRJNA722100_discovery_only"
conflicts <- weights[weights$weight == 0, , drop = FALSE]
if (nrow(conflicts)) {
  conflict_entries <- merge(
    signature[, c("gene", "cell_type", "rank_within_lineage",
                  "mean_donor_effect", "weight")],
    conflicts[, c("gene", "n_lineages")], by = "gene"
  )
} else {
  conflict_entries <- data.frame(
    gene = character(), cell_type = character(), rank_within_lineage = integer(),
    mean_donor_effect = numeric(), weight = numeric(), n_lineages = integer()
  )
}
write_tsv(conflict_entries,
          "results/tables/bulk/frozen_multicellular_signature_direction_conflicts.tsv")
weights <- weights[weights$weight != 0, , drop = FALSE]
stopifnot(nrow(weights) <= nrow(signature), nrow(weights) >= 140L,
          !any(weights$weight == 0))

write_tsv(signature, "data/processed/frozen_multicellular_signature.tsv")
write_tsv(weights, "data/processed/frozen_multicellular_signature_gene_weights.tsv")
write_tsv(signature, "results/tables/bulk/frozen_multicellular_signature.tsv")
write_tsv(data.frame(
  discovery_dataset = "PRJNA722100",
  n_lineages = length(unique(signature$cell_type)),
  lineage_gene_entries = nrow(signature),
  unique_genes_before_direction_conflict_exclusion = length(unique(signature$gene)),
  direction_conflict_genes_excluded = nrow(conflicts),
  usable_unique_genes = nrow(weights),
  validation_expression_used_for_selection = FALSE,
  validation_labels_used_for_selection = FALSE
), "results/tables/bulk/frozen_multicellular_signature_manifest.tsv")

message("Fixed discovery-derived signature: 150 lineage-gene entries, ",
        nrow(weights), " usable unique genes; ", nrow(conflicts),
        " direction-conflict gene(s) excluded.")
