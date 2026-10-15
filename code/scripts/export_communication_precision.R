#!/usr/bin/env Rscript
# Optional read-only source export from the original four native LIANA caches.
# No inference is run. Only serialization-boundary values and source hashes are
# released; the original cell objects/native caches remain outside the package.
args <- commandArgs(TRUE)
cache_arg <- grep("^--cache-dir=", args, value = TRUE)
if (length(cache_arg) != 1L) stop("Supply --cache-dir=/path/to/original/communication_cache")
cache_dir <- sub("^--cache-dir=", "", cache_arg)
output_arg <- grep("^--output-dir=", args, value = TRUE)
out <- if (length(output_arg)) sub("^--output-dir=", "", output_arg[[1]]) else "analysis"
sha <- function(path) {
  command <- if (nzchar(Sys.which("sha256sum"))) "sha256sum" else "shasum"
  options <- if (command == "shasum") c("-a", "256") else character()
  value <- system2(command, c(options, shQuote(path)), stdout = TRUE)
  if (!is.null(attr(value, "status"))) stop("Cannot hash source: ", path)
  substr(value[[1]], 1, 64)
}
published_path <- "results/tables/communication/LIANA_NATMI_CellPhoneDB_all.tsv.gz"
published <- read.delim(gzfile(published_path), check.names = FALSE, stringsAsFactors = FALSE)
samples <- c("H01045__concave", "H01045__convex", "H66206__concave", "H66206__convex")
all_native <- list(); manifests <- list()
for (sample in samples) {
  filename <- paste0(sample, "__soupx_consensus_v2_20260731.rds")
  path <- file.path(cache_dir, filename)
  native <- readRDS(path)
  stopifnot(identical(sort(names(native)), sort(paste(sample, c("natmi", "cellphonedb")))))
  for (name in names(native)) {
    d <- native[[name]]
    stopifnot(is.data.frame(d), all(d$sample_id == sample))
    d$cache_basename <- filename; d$cache_sha256 <- sha(path)
    all_native[[name]] <- d
  }
  manifests[[sample]] <- data.frame(sample_id = sample, cache_basename = filename,
    cache_bytes = file.info(path)$size, cache_sha256 = sha(path),
    cache_shape = "Legacy list of two native method data.frames; no embedded provenance stamp",
    export_scope = "Read-only native-cache comparison and precision export; no LIANA rerun or retrospective original input binding")
}
native <- lapply(all_native, function(d) {
  for (column in setdiff(names(published), names(d))) d[[column]] <- NA
  d[, c(names(published), "cache_basename", "cache_sha256")]
})
native <- do.call(rbind, native)
keys <- c("method", "sample_id", "source", "target", "ligand.complex", "ligand", "receptor.complex", "receptor")
key <- function(d) do.call(paste, c(d[keys], sep = "\034"))
a <- key(native); b <- key(published)
stopifnot(nrow(native) == nrow(published), nrow(native) == 144872L,
          !anyDuplicated(a), !anyDuplicated(b), setequal(a, b))
native <- native[match(b, a), , drop = FALSE]
max_error <- 0
for (column in names(published)) {
  x <- native[[column]]; y <- published[[column]]
  stopifnot(identical(is.na(x), is.na(y)))
  ok <- !is.na(x)
  if (is.numeric(x) && is.numeric(y)) {
    error <- abs(x[ok] - y[ok]); max_error <- max(max_error, error)
    if (any(error > 1e-12 + 1e-10 * abs(y[ok]))) stop("Native/published numerical mismatch: ", column)
  } else if (!identical(as.character(x[ok]), as.character(y[ok]))) {
    stop("Native/published label/flag mismatch: ", column)
  }
}
support <- ifelse(native$method == "natmi", TRUE, native$pvalue <= .05)
stopifnot(identical(support, native$method_supported))
boundary <- published$method == "cellphonedb" & !is.na(published$pvalue) & published$pvalue == .05
stopifnot(sum(boundary) == 37L, all(!support[boundary]), all(native$pvalue[boundary] > .05))
result <- native[boundary, c(keys, "patient_id", "side", "cache_basename", "cache_sha256")]
result$published_row <- which(boundary)
result$published_pvalue_text <- "0.05"
result$native_pvalue_17g <- sprintf("%.17g", native$pvalue[boundary])
result$method_supported <- support[boundary]
result$published_all_table_sha256 <- sha(published_path)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.table(result, file.path(out, "communication_native_precision_boundaries.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
manifest <- do.call(rbind, manifests)
manifest$published_all_table_sha256 <- sha(published_path)
manifest$published_rows_matched <- nrow(published)
manifest$maximum_native_vs_printed_absolute_error <- sprintf("%.17g", max_error)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
stopifnot(length(script_arg) == 1L)
manifest$exporter_sha256 <- sha(sub("^--file=", "", script_arg))
manifest$comparison_keys <- paste(keys, collapse = ";")
write.table(manifest, file.path(out, "communication_native_precision_sources.tsv"), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")
cat("All native/published method records and flags matched:", nrow(published),
    "rows. Exported", nrow(result), "actual native boundary P values.\n")
