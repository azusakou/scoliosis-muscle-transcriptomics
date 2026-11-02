#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(scuttle); library(liana); library(Matrix)
})
source("scripts/lib_common.R")
set.seed(501422)

obj <- readRDS("data/processed/scrna/scrna_annotated.rds")
obj <- JoinLayers(obj, assay = "SoupX")
counts <- LayerData(obj, assay = "SoupX", layer = "counts")
meta <- obj@meta.data
all_results <- list()
CACHE_VERSION <- "soupx_consensus_v3_input_checked"
validate_paired_samples(unique(meta[, c("sample_id", "patient_id", "side")]), expected_pairs = 2)
stopifnot(!anyNA(meta$cell_type))
cache_provenance <- list(
  input_md5 = unname(tools::md5sum("data/processed/scrna/scrna_annotated.rds")),
  script_md5 = unname(tools::md5sum("scripts/50_cell_communication.R")),
  versions = vapply(c("liana", "scuttle", "SingleCellExperiment", "Matrix"),
                    function(p) as.character(packageVersion(p)), character(1)),
  liana_source_commit = utils::packageDescription("liana")$RemoteSha,
  methods = c("natmi", "cellphonedb"), resource = "Consensus", min_cells = 20,
  count_assay = "SoupX")
dir.create("data/processed/scrna/communication_cache", recursive = TRUE, showWarnings = FALSE)
for (sid in unique(meta$sample_id)) {
  cache_file <- file.path(
    "data/processed/scrna/communication_cache",
    paste0(sid, "__", CACHE_VERSION, ".rds")
  )
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (!identical(cached$provenance, cache_provenance) || !identical(cached$sample_id, sid)) {
      stop("Communication cache input/version mismatch: ", cache_file,
           ". Preserve/move this cache and recompute.")
    }
    sample_results <- cached$results
    all_results <- c(all_results, sample_results)
    next
  }
  ix <- meta$sample_id == sid
  sce <- SingleCellExperiment(assays = list(counts = counts[, ix, drop = FALSE]), colData = meta[ix, , drop = FALSE])
  sce <- logNormCounts(sce)
  result <- liana_wrap(sce, method = c("natmi", "cellphonedb"), resource = "Consensus",
                       idents_col = "cell_type", min_cells = 20, verbose = FALSE)
  sample_results <- list()
  for (method in names(result)) {
    z <- as.data.frame(result[[method]])
    z$method <- method; z$sample_id <- sid
    z$patient_id <- unique(meta$patient_id[ix]); z$side <- unique(meta$side[ix])
    z$count_assay <- "SoupX_ambient_corrected_counts"
    z$cache_version <- CACHE_VERSION
    z$method_score <- if (method == "natmi") z$edge_specificity else z$lr.mean
    z$method_supported <- if (method == "natmi") TRUE else z$pvalue <= 0.05
    # Preserve only this probability at round-trip precision in future raw runs.
    # The original published all-method gzip is immutable; its actual native
    # boundary values are separately source-bound in the retained release.
    z$native_pvalue_17g <- if (method == "natmi") NA_character_ else sprintf("%.17g", z$pvalue)
    sample_results[[paste(sid, method)]] <- z
  }
  temporary <- tempfile(tmpdir = dirname(cache_file))
  saveRDS(list(provenance = cache_provenance, sample_id = sid, results = sample_results), temporary)
  if (!file.rename(temporary, cache_file)) stop("Cannot install completed cache: ", cache_file)
  all_results <- c(all_results, sample_results)
}
all_results <- dplyr::bind_rows(all_results)
write_tsv(all_results, "results/tables/communication/LIANA_NATMI_CellPhoneDB_all.tsv.gz")

# Require detectable ligand/receptor, cross-lineage edges, and paired donor direction in each method.
z <- all_results[all_results$ligand.prop >= 0.1 & all_results$receptor.prop >= 0.1 &
                   all_results$source != all_results$target & all_results$method_supported, ]
key <- c("method", "source", "target", "ligand", "receptor", "patient_id")
conc <- z[z$side == "concave", c(key, "method_score")]
conv <- z[z$side == "convex", c(key, "method_score")]
paired <- merge(conc, conv, by = key, suffixes = c(".concave", ".convex"))
paired$paired_score_difference <- paired$method_score.concave - paired$method_score.convex
summary_method <- aggregate(paired_score_difference ~ method + source + target + ligand + receptor, paired,
                            function(x) c(mean = mean(x), donor_concordant = as.numeric(length(x)==2 && (all(x>0)||all(x<0))), n=length(x)))
summary_method <- cbind(summary_method[1:5], as.data.frame(summary_method$paired_score_difference))

nat <- summary_method[summary_method$method == "natmi", ]; cp <- summary_method[summary_method$method == "cellphonedb", ]
cross <- merge(nat, cp, by = c("source","target","ligand","receptor"), suffixes = c(".natmi", ".cellphonedb"))
cross$method_direction_concordant <- sign(cross$mean.natmi) == sign(cross$mean.cellphonedb)

# Check the same direction against ligand+receptor donor-level pseudobulk effects.
pb <- read.delim(gzfile("results/tables/pseudobulk/scrna_donor_specific_effects.tsv.gz"), check.names = FALSE)
lookup <- function(gene, cell_type, patient) {
  q <- pb$paired_logCPM_difference[pb$gene == gene & pb$cell_type == cell_type & pb$patient_id == patient]
  if (length(q)) q[1] else NA_real_
}
for (pid in unique(pb$patient_id)) {
  cross[[paste0("pseudobulk_joint_effect_", pid)]] <- mapply(function(lig, rec, src, tgt)
    lookup(lig, src, pid) + lookup(rec, tgt, pid), cross$ligand, cross$receptor, cross$source, cross$target)
}
joint_cols <- grep("^pseudobulk_joint_effect", names(cross), value = TRUE)
cross$pseudobulk_donor_concordant <- apply(cross[, joint_cols, drop = FALSE], 1,
                                           function(x) all(is.finite(x)) && (all(x>0)||all(x<0)))
cross$pseudobulk_matches_methods <- cross$pseudobulk_donor_concordant &
  sign(rowMeans(cross[, joint_cols, drop = FALSE], na.rm = TRUE)) == sign(cross$mean.natmi)
cross$retained <- cross$donor_concordant.natmi == 1 & cross$donor_concordant.cellphonedb == 1 &
  cross$method_direction_concordant & cross$pseudobulk_matches_methods
cross$inferential_status <- "exploratory_n2_requires_external_validation"
cross <- cross[order(cross$retained, abs(cross$mean.natmi) + abs(cross$mean.cellphonedb), decreasing = TRUE), ]
write_tsv(paired, "results/tables/communication/paired_method_edges_by_donor.tsv.gz")
write_tsv(cross, "results/tables/communication/cross_method_donor_aware_edges.tsv")

capture_versions(c("liana", "SingleCellExperiment", "scuttle", "Seurat"), "logs/communication_package_versions.tsv")
message("SoupX communication complete: ", sum(cross$retained),
        " edges passed both methods and donor-aware filters.")
