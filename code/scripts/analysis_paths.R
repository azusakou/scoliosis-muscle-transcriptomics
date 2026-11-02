# Read-only aliases for byte-identical published analysis inputs.  A raw rebuild
# writes its own legacy output path; that fresh output takes precedence without
# copying or overwriting the separately retained, checksum-bound baseline.
analysis_input_aliases <- c(
  "results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz" = "analysis/bulk_gate_inputs/original_A_gene_results.tsv.gz",
  "results/tables/bulk_gse311180_requant/GSE311180_paired_edgeR_full.tsv.gz" = "analysis/bulk_gate_inputs/original_B_gene_results.tsv.gz",
  "results/tables/bulk_gse311180_requant/GSE311180_QC.tsv" = "analysis/bulk_gate_inputs/original_B_QC.tsv",
  "analysis/bulk_edgeR_null_logFC_GSE254300.tsv.gz" = "analysis/bulk_gate_inputs/A32_logFC.tsv.gz",
  "analysis/bulk_edgeR_null_grid_GSE254300.tsv" = "analysis/bulk_gate_inputs/A32_grid.tsv",
  "analysis/audit_gse311180_rna_expected_pairs.tsv" = "analysis/bulk_gate_inputs/RNA_expected_pairs.tsv",
  "analysis/bulk_gate_normalization_factors.tsv" = "analysis/bulk_gate_inputs/published_retained_TMM_factors.tsv",
  "results/figures/source_data/Figure_S3_osteoblast_control.tsv" = "results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_paired_effects.tsv"
)

analysis_input <- function(path, root = ".") {
  if (length(path) != 1L || is.na(path) || !nzchar(path) ||
      grepl("^[/~]|^[A-Za-z]:|\\\\", path) ||
      any(strsplit(path, "/", fixed = TRUE)[[1]] %in% c("", ".", "..")))
    stop("Unsafe analysis input path")
  root <- normalizePath(root, mustWork = TRUE)
  checked <- function(relative) {
    parts <- strsplit(relative, "/", fixed = TRUE)[[1]]
    candidate <- root
    for (part in parts) {
      candidate <- file.path(candidate, part)
      link <- Sys.readlink(candidate)
      if (!is.na(link) && nzchar(link)) stop("Symlink analysis input path: ", relative)
    }
    if (!file.exists(candidate)) return(NULL)
    actual <- normalizePath(candidate, mustWork = TRUE)
    if (dir.exists(actual) || !startsWith(actual, paste0(root, .Platform$file.sep)))
      stop("Analysis input escapes root or is not a regular file")
    actual
  }
  fresh <- checked(path)
  canonical <- unname(analysis_input_aliases[path])
  if (!is.na(canonical)) {
    # Validate even a non-selected baseline: a symlink is never a valid alias.
    fixed <- checked(canonical)
    if (!is.null(fresh)) {
      message("Analysis input scope: fresh producer output ", path,
              "; immutable baseline ", canonical, " is not substituted.")
      return(fresh)
    }
    if (!is.null(fixed)) return(fixed)
  } else if (!is.null(fresh)) return(fresh)
  stop("Missing analysis input: ", path)
}
