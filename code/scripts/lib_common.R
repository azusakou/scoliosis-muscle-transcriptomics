suppressPackageStartupMessages({
  library(ggplot2)
})

source("scripts/analysis_paths.R")

raw_data_root <- function() Sys.getenv("AIS_RAW_ROOT", "data/raw/public")

validate_annotation_map <- function(mapping, observed_clusters) {
  required <- c("cluster", "cell_type", "cell_state_label", "canonical_marker_evidence", "confidence")
  if (!all(required %in% names(mapping)) || !nrow(mapping) || anyNA(mapping[, required]))
    stop("Missing or incomplete annotation-map fields")
  clusters <- as.character(mapping$cluster)
  if (anyDuplicated(clusters) || any(!nzchar(trimws(clusters))))
    stop("Duplicate or empty annotation-map cluster IDs")
  if (any(vapply(mapping[, required[-1], drop = FALSE],
                 function(x) any(!nzchar(trimws(as.character(x)))), logical(1))))
    stop("Empty annotation label, marker evidence or confidence")
  if (!all(mapping$confidence %in% c("high", "medium", "low")))
    stop("Unexpected annotation confidence category")
  observed <- unique(as.character(observed_clusters))
  if (!length(observed) || anyNA(observed) || any(!nzchar(trimws(observed))) ||
      !setequal(observed, clusters) || length(observed) != nrow(mapping))
    stop("Annotation map does not cover each observed cluster exactly once")
  invisible(TRUE)
}

validate_paired_samples <- function(info, sample_column = "sample_id", expected_pairs = NULL) {
  required <- c(sample_column, "patient_id", "side")
  if (!all(required %in% names(info)) || !nrow(info) || anyNA(info[, required]))
    stop("Missing or incomplete patient/sample/side metadata")
  if (anyDuplicated(info[[sample_column]]) || any(info[[sample_column]] == "") ||
      any(info$patient_id == "") || any(!info$side %in% c("concave", "convex")))
    stop("Duplicate sample IDs, empty patient IDs or unexpected side labels")
  pairs <- table(info$patient_id, factor(info$side, levels = c("concave", "convex")))
  if (any(pairs != 1L)) stop("Each patient must have exactly one concave and one convex library")
  if (!is.null(expected_pairs) && nrow(pairs) != expected_pairs)
    stop("Expected ", expected_pairs, " patient pairs; found ", nrow(pairs))
  invisible(TRUE)
}


composition_from_cells <- function(metadata) {
  required <- c("sample_id", "patient_id", "side", "cell_type")
  if (!all(required %in% names(metadata)) || !nrow(metadata) || anyNA(metadata[, required]) ||
      any(metadata$cell_type == "")) stop("Invalid cell composition metadata")
  samples <- unique(metadata[, c("sample_id", "patient_id", "side")])
  validate_paired_samples(samples)
  counts <- as.data.frame(table(sample_id = metadata$sample_id, cell_type = metadata$cell_type),
                          stringsAsFactors = FALSE)
  # Every observed lineage has a row in every library, including true zero cells.
  counts <- merge(counts, samples, by = "sample_id", sort = TRUE)
  totals <- aggregate(Freq ~ sample_id, counts, sum)
  names(totals)[2] <- "sample_total"
  counts <- merge(counts, totals, by = "sample_id", sort = TRUE)
  counts$proportion <- counts$Freq / counts$sample_total
  counts[, c("sample_id", "patient_id", "side", "cell_type", "Freq", "sample_total", "proportion")]
}

palette_contract <- c(
  concave = "#D24B40",
  convex = "#3182BD",
  neutral_dark = "#272727",
  neutral_mid = "#767676",
  neutral_light = "#D8D8D8",
  signal_teal = "#33B5A5",
  accent_orange = "#E28E2C"
)

theme_nature_contract <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.7),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      panel.grid = element_blank()
    )
}

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 120, dpi = 600) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite::svglite(paste0(filename, ".svg"), width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()
  grDevices::pdf(paste0(filename, ".pdf"), width = width_in, height = height_in,
                 family = "Helvetica", useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
  ragg::agg_tiff(paste0(filename, ".tiff"), width = width_in, height = height_in,
                 units = "in", res = dpi, compression = "lzw")
  print(plot)
  grDevices::dev.off()
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else path
  on.exit(if (inherits(connection, "connection")) close(connection), add = TRUE)
  # Plot labels can contain newlines. Quote only columns requiring escaping,
  # preserving the historical representation of ordinary numerical tables.
  columns <- if (is.data.frame(x)) x else as.data.frame(x)
  unsafe <- vapply(columns, function(v) any(grepl('[\t\r\n"]', as.character(v)), na.rm = TRUE), logical(1))
  utils::write.table(x, connection, sep = "\t", quote = which(unsafe), qmethod = "double",
                     row.names = FALSE, na = "NA")
}

robust_qc_thresholds <- function(total, detected, percent_mt) {
  log_total <- log10(pmax(total, 1))
  log_detected <- log10(pmax(detected, 1))
  safe_mad <- function(x) {
    value <- stats::mad(x, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(value) || value == 0) value <- stats::IQR(x, na.rm = TRUE) / 1.349
    if (!is.finite(value) || value == 0) value <- 0.1
    value
  }
  list(
    min_counts = floor(max(100, 10^(stats::median(log_total) - 3 * safe_mad(log_total)))),
    max_counts = ceiling(10^(stats::median(log_total) + 5 * safe_mad(log_total))),
    min_features = floor(max(200, 10^(stats::median(log_detected) - 3 * safe_mad(log_detected)))),
    max_features = ceiling(10^(stats::median(log_detected) + 5 * safe_mad(log_detected))),
    max_percent_mt = min(40, stats::median(percent_mt, na.rm = TRUE) + 3 * safe_mad(percent_mt))
  )
}

capture_versions <- function(packages, path) {
  versions <- vapply(packages, function(package) {
    if (package == "R") as.character(getRversion()) else if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else "MISSING"
  }, character(1))
  write_tsv(data.frame(package = names(versions), version = unname(versions)), path)
}

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  entries <- strsplit(lines, "\t", fixed = FALSE)
  stats::setNames(lapply(entries, function(x) unique(x[-c(1, 2)])),
                  vapply(entries, `[[`, character(1), 1))
}

gene_set_zscore <- function(log_expression, gene_sets, min_genes = 5) {
  gene_z <- t(scale(t(log_expression)))
  gene_z[!is.finite(gene_z)] <- 0
  used <- lapply(gene_sets, intersect, y = rownames(gene_z))
  used <- used[lengths(used) >= min_genes]
  if (!length(used)) return(matrix(numeric(), nrow = 0, ncol = ncol(log_expression)))
  scores <- vapply(used, function(g) colMeans(gene_z[g, , drop = FALSE]), numeric(ncol(gene_z)))
  scores <- t(scores)
  colnames(scores) <- colnames(log_expression)
  scores
}

# Build a projection matrix without using condition labels. Differential-expression
# filters may depend on the contrast and must therefore never define a validation
# signature. Here detectability depends only on CPM and the number of libraries.
label_blind_logcpm <- function(counts, norm_factors, min_cpm = 1,
                               min_samples = max(2L, ceiling(0.2 * ncol(counts)))) {
  stopifnot(length(norm_factors) == ncol(counts), min_samples <= ncol(counts))
  if (any(!is.finite(norm_factors)) || any(norm_factors <= 0) ||
      length(min_cpm) != 1L || !is.finite(min_cpm) || min_cpm <= 0 ||
      length(min_samples) != 1L || min_samples < 1L || min_samples != as.integer(min_samples))
    stop("Invalid normalization factors or detectability threshold")
  y <- edgeR::DGEList(counts)
  y$samples$norm.factors <- norm_factors
  raw_cpm <- edgeR::cpm(y, log = FALSE)
  detected <- rowSums(raw_cpm >= min_cpm) >= min_samples
  logcpm <- edgeR::cpm(y[detected, , keep.lib.sizes = TRUE],
                       log = TRUE, prior.count = 2)
  list(logcpm = logcpm, detected = detected, min_cpm = min_cpm,
       min_samples = min_samples)
}

project_frozen_signature <- function(log_expression, signature, min_genes = 20L) {
  if (is.null(rownames(log_expression)) || anyDuplicated(rownames(log_expression)) ||
      is.null(colnames(log_expression)) || anyDuplicated(colnames(log_expression)) ||
      ncol(log_expression) < 2L || any(!is.finite(log_expression)))
    stop("Projection requires finite expression with unique gene and library identifiers")
  if (!all(c("gene", "weight") %in% names(signature)) || anyNA(signature[, c("gene", "weight")]) ||
      any(signature$gene == "") || any(!is.finite(signature$weight)))
    stop("Missing or invalid signature identities/weights")
  if ("cell_type" %in% names(signature) && anyDuplicated(signature[, c("gene", "cell_type")]))
    stop("Duplicated gene-lineage signature entry")
  if (length(min_genes) != 1L || !is.finite(min_genes) || min_genes < 1L || min_genes != as.integer(min_genes))
    stop("min_genes must be a positive integer")
  weights <- stats::aggregate(weight ~ gene, signature, mean)
  weights <- weights[is.finite(weights$weight) & weights$weight != 0, , drop = FALSE]
  weights$available <- weights$gene %in% rownames(log_expression)
  used <- weights[weights$available, , drop = FALSE]
  if (nrow(used) < min_genes) stop("Fewer than ", min_genes,
                                  " frozen genes passed label-blind detectability")
  z <- t(scale(t(log_expression[used$gene, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  score <- colSums(z * used$weight) / sum(abs(used$weight))
  list(score = score, weights = weights, used_weights = used, gene_z = z)
}

sample_matched_genes_without_replacement <- function(target_genes, target_bins,
                                                     universe_bins,
                                                     excluded_genes = target_genes) {
  stopifnot(length(target_genes) == length(target_bins),
            !is.null(names(universe_bins)))
  if (anyNA(target_bins) || anyNA(universe_bins) || anyDuplicated(names(universe_bins)) ||
      anyDuplicated(target_genes)) stop("Invalid matching bins or duplicate gene identities")
  sampled <- character(length(target_genes))
  for (bin in unique(target_bins)) {
    positions <- which(target_bins == bin)
    pool <- names(universe_bins)[universe_bins == bin &
                                  !names(universe_bins) %in% excluded_genes]
    if (length(pool) < length(positions)) {
      stop("Insufficient unique genes in expression bin ", bin,
           ": need ", length(positions), ", found ", length(pool))
    }
    sampled[positions] <- sample(pool, length(positions), replace = FALSE)
  }
  stopifnot(length(unique(sampled)) == length(sampled))
  sampled
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
