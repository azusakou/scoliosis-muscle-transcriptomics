#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(edgeR); library(fgsea) })
source("scripts/lib_common.R")

pb <- readRDS("data/processed/scrna/pseudobulk_logcpm_by_cell_type.rds")
hallmark <- read_gmt("references/gene_sets/h.all.v2026.1.Hs.symbols.gmt")
reactome <- read_gmt("references/gene_sets/c2.cp.reactome.v2026.1.Hs.symbols.gmt")
tft <- read_gmt("references/gene_sets/c3.tft.gtrd.v2026.1.Hs.symbols.gmt")

hallmark_keep <- c(
  "HALLMARK_MYOGENESIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_ANGIOGENESIS",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_ESTROGEN_RESPONSE_EARLY",
  "HALLMARK_ESTROGEN_RESPONSE_LATE", "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_NOTCH_SIGNALING", "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_MTORC1_SIGNALING"
)
reactome_patterns <- c("EXTRACELLULAR_MATRIX_ORGANIZATION", "COLLAGEN_FORMATION", "INTEGRIN_CELL_SURFACE_INTERACTIONS",
                       "SIGNALING_BY_NOTCH", "SIGNALING_BY_WNT", "SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",
                       "MITOCHONDRIAL_BIOGENESIS", "RESPIRATORY_ELECTRON_TRANSPORT", "VEGF", "JAK_STAT")
reactome_keep <- unique(unlist(lapply(reactome_patterns, function(p) grep(p, names(reactome), value = TRUE))))
pathways <- c(hallmark[intersect(hallmark_keep, names(hallmark))], reactome[reactome_keep])

score_long <- list(); effect_long <- list(); tf_long <- list(); overlap <- list()
for (ct in names(pb)) {
  mat <- pb[[ct]]$logcpm
  info <- pb[[ct]]$sample_info
  scores <- gene_set_zscore(mat, pathways, min_genes = 5)
  tf_scores <- gene_set_zscore(mat, tft, min_genes = 10)
  if (nrow(tf_scores) > 0) {
    v <- apply(tf_scores, 1, var)
    tf_scores <- tf_scores[order(v, decreasing = TRUE)[seq_len(min(100, length(v)))], , drop = FALSE]
  }
  score_long[[ct]] <- do.call(rbind, lapply(rownames(scores), function(gs) data.frame(
    cell_type = ct, gene_set = gs, pb_id = colnames(scores), score = scores[gs, ],
    patient_id = info$patient_id[match(colnames(scores), info$pb_id)],
    side = as.character(info$side[match(colnames(scores), info$pb_id)]))))
  if (nrow(tf_scores)) tf_long[[ct]] <- do.call(rbind, lapply(rownames(tf_scores), function(gs) data.frame(
    cell_type = ct, tf_target_set = gs, pb_id = colnames(tf_scores), score = tf_scores[gs, ],
    patient_id = info$patient_id[match(colnames(tf_scores), info$pb_id)],
    side = as.character(info$side[match(colnames(tf_scores), info$pb_id)]))))
  overlap[[ct]] <- data.frame(cell_type = ct, gene_set = names(pathways),
                              n_set = lengths(pathways), n_overlap = lengths(lapply(pathways, intersect, rownames(mat))))
}
scores <- do.call(rbind, score_long)
effects <- merge(scores[scores$side == "concave", ], scores[scores$side == "convex", ],
                 by = c("cell_type", "gene_set", "patient_id"), suffixes = c(".concave", ".convex"))
effects$paired_difference <- effects$score.concave - effects$score.convex
summary_effects <- aggregate(paired_difference ~ cell_type + gene_set, effects,
                             function(x) c(mean = mean(x), concordant = as.numeric(all(x>0)||all(x<0))))
summary_effects <- cbind(summary_effects[c("cell_type", "gene_set")], as.data.frame(summary_effects$paired_difference))
donor_wide <- reshape(
  effects[, c("cell_type", "gene_set", "patient_id", "paired_difference")],
  idvar = c("cell_type", "gene_set"), timevar = "patient_id", direction = "wide"
)
summary_effects <- merge(summary_effects, donor_wide, by = c("cell_type", "gene_set"), all.x = TRUE)
summary_effects$inferential_status <- "descriptive_n2"
write_tsv(scores, "results/tables/programs/scrna_pathway_scores.tsv")
write_tsv(effects, "results/tables/programs/scrna_pathway_paired_effects_by_donor.tsv")
write_tsv(summary_effects, "results/tables/programs/scrna_pathway_effect_summary.tsv")
write_tsv(do.call(rbind, tf_long), "results/tables/programs/scrna_gtrd_tf_target_scores_top_variable.tsv")
write_tsv(do.call(rbind, overlap), "results/tables/programs/gene_set_overlap_scrna.tsv")
saveRDS(list(pathways = pathways, hallmark_release = "2026.1.Hs", scoring = "mean_per_gene_zscore"),
        "data/processed/frozen_pathway_definitions.rds")
message("Frozen pathway scores complete.")
