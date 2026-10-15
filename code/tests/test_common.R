#!/usr/bin/env Rscript
source("scripts/lib_common.R")
checked <- character()
check <- function(label, expression) {
  force(expression)
  checked <<- c(checked, label)
  cat("PASS:", label, "\n")
}
errors <- function(expression) inherits(tryCatch({force(expression); NULL}, error = identity), "error")
check("annotation map has one complete assignment per observed cluster", {
  mapping <- read.delim("config/annotation_cluster_map.tsv", check.names = FALSE)
  validate_annotation_map(mapping, as.character(0:16))
  stopifnot(errors(validate_annotation_map(rbind(mapping, mapping[1, ]), as.character(0:16))),
    errors(validate_annotation_map(mapping[-1, ], as.character(0:16))))
  bad <- mapping; bad$cell_type[1] <- ""
  stopifnot(errors(validate_annotation_map(bad, as.character(0:16))))
  bad <- mapping; bad$confidence[1] <- NA
  stopifnot(errors(validate_annotation_map(bad, as.character(0:16))))
  bad <- mapping; bad$confidence[1] <- "unspecified"
  stopifnot(errors(validate_annotation_map(bad, as.character(0:16))))
})
check("one library on each side per patient", {
  metadata <- data.frame(sample_id = c("p1a", "p1b", "p2a", "p2b"),
    patient_id = rep(c("p1", "p2"), each = 2), side = rep(c("concave", "convex"), 2))
  validate_paired_samples(metadata, expected_pairs = 2)
  stopifnot(errors(validate_paired_samples(metadata[-1, ])),
    errors(validate_paired_samples(rbind(metadata, metadata[1, ]))),
    errors(validate_paired_samples(transform(metadata, side = "concave"))))
})
check("zero-cell lineages remain present in each library", {
  cells <- data.frame(sample_id = c("p1a", "p1a", "p1b"), patient_id = "p1",
                      side = c("concave", "concave", "convex"), cell_type = c("a", "b", "a"))
  composition <- composition_from_cells(cells)
  zero <- composition[composition$sample_id == "p1b" & composition$cell_type == "b", ]
  stopifnot(nrow(composition) == 4, nrow(zero) == 1, zero$Freq == 0, zero$proportion == 0)
})
check("signed projection, conflict removal and constant rows", {
  expression <- rbind(g1 = c(1, 2, 3, 4), g2 = c(4, 3, 2, 1), flat = c(9, 9, 9, 9),
                      conflict = c(1, 2, 3, 4))
  colnames(expression) <- paste0("s", 1:4)
  signature <- data.frame(gene = c("g1", "g2", "flat", "conflict", "conflict"),
    weight = c(1, -1, 1, -1, 1), cell_type = c("a", "a", "a", "a", "b"))
  result <- project_frozen_signature(expression, signature, min_genes = 3)
  expected <- (expression["g1", ] - 2.5) / sqrt(5 / 3) * 2 / 3
  stopifnot(max(abs(result$score - expected)) < 1e-14,
    nrow(result$used_weights) == 3, all(result$gene_z["flat", ] == 0),
    !"conflict" %in% result$weights$gene,
    errors(project_frozen_signature(expression, signature, min_genes = 4)))
  bad <- expression; rownames(bad)[1] <- rownames(bad)[2]
  stopifnot(errors(project_frozen_signature(bad, signature, min_genes = 1)))
  bad <- expression; bad[1, 1] <- NA_real_
  stopifnot(errors(project_frozen_signature(bad, signature, min_genes = 1)),
    errors(project_frozen_signature(expression, rbind(signature, signature[1, ]), min_genes = 1)))
})
check("label-blind detectability and retained library sizes", {
  counts <- rbind(common = c(1000, 1100, 1200, 1300), absent = c(0, 0, 0, 0),
                  rare = c(10, 0, 0, 0))
  colnames(counts) <- paste0("s", 1:4)
  nf <- rep(1, 4)
  result <- label_blind_logcpm(counts, nf, min_cpm = 1, min_samples = 2)
  y <- edgeR::DGEList(counts); y$samples$norm.factors <- nf
  expected <- edgeR::cpm(y["common", , keep.lib.sizes = TRUE], log = TRUE, prior.count = 2)
  stopifnot(identical(as.logical(result$detected), c(TRUE, FALSE, FALSE)),
    max(abs(result$logcpm - expected)) < 1e-14,
    errors(label_blind_logcpm(counts, c(1, -1, 1, 1))),
    errors(label_blind_logcpm(counts, nf, min_samples = 0)),
    errors(label_blind_logcpm(counts, nf, min_samples = 1.5)))
})
check("matched sampling without replacement or signature contamination", {
  bins <- setNames(rep(1:2, each = 10), paste0("g", 1:20))
  targets <- c("g1", "g2", "g11", "g12")
  set.seed(1)
  for (i in 1:100) {
    sampled <- sample_matched_genes_without_replacement(targets, bins[targets], bins)
    stopifnot(!anyDuplicated(sampled), !any(sampled %in% targets), identical(unname(bins[sampled]), unname(bins[targets])))
  }
  stopifnot(errors(sample_matched_genes_without_replacement(targets, bins[targets], bins[targets])))
})
check("TSV tab, quote and multiline labels round-trip", {
  directory <- tempfile(); dir.create(directory)
  original <- data.frame(label = c("two\nlines", "a\tb", 'say "hello"'), value = c(.1, .2, .3))
  for (extension in c(".tsv", ".tsv.gz")) {
    path <- file.path(directory, paste0("table", extension)); write_tsv(original, path)
    connection <- if (extension == ".tsv.gz") gzfile(path, "rt") else file(path, "rt")
    restored <- read.delim(connection, check.names = FALSE); close(connection)
    stopifnot(identical(original, restored))
  }
  unlink(directory, recursive = TRUE)
})
check("synthetic per-cell before/after audit uses actual counts and cells", {
  suppressPackageStartupMessages(library(Seurat))
  directory <- tempfile(); dir.create(directory)
  cells <- paste0("cell", seq_len(14846))
  counts <- Matrix::Matrix(rbind(IGKC = rep(c(0, 2, 3), length.out = 14846), ACTA1 = rep(5, 14846)), sparse = TRUE)
  colnames(counts) <- cells
  object <- CreateSeuratObject(counts = counts)
  corrected <- counts; corrected["IGKC", ] <- 0
  object[["SoupX"]] <- CreateAssay5Object(counts = corrected)
  object$sample_id <- rep(c("p1__concave", "p1__convex", "p2__concave", "p2__convex"), length.out = 14846)
  object$patient_id <- sub("__.*", "", object$sample_id)
  object$side <- sub(".*__", "", object$sample_id)
  object$cell_type <- rep(paste0("lineage", 1:9), length.out = 14846)
  input <- file.path(directory, "synthetic.rds"); saveRDS(object, input)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("analysis/signature_cell_expression_audit.R", shQuote(input), shQuote(directory)),
    stdout = file.path(directory, "run.log"), stderr = file.path(directory, "run.log"))
  if (status != 0) stop(paste(readLines(file.path(directory, "run.log")), collapse = "\n"))
  table <- read.delim(file.path(directory, "signature_expression_by_library_lineage_before_after.tsv"))
  stopifnot(nrow(table) == 144)
  for (i in seq_len(nrow(table))) {
    row <- table[i, ]; ix <- which(object$sample_id == row$sample_id & object$cell_type == row$cell_type)
    values <- if (row$assay == "RNA") counts[row$gene, ix] else corrected[row$gene, ix]
    stopifnot(row$n_cells == length(ix), row$n_detected == sum(values > 0),
      abs(row$detection_fraction - mean(values > 0)) < 1e-12,
      abs(row$mean_UMI_per_cell - mean(values)) < 1e-12)
  }
  availability <- read.delim(file.path(directory, "signature_requested_gene_availability.tsv"))
  stopifnot(all(!availability$available[availability$gene == "ESR1"]))
  unlink(directory, recursive = TRUE)
})
cat("All", length(checked), "R regression groups passed. Synthetic data are tests only.\n")
