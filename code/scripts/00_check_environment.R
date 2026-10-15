#!/usr/bin/env Rscript

required <- c(
  "BiocParallel", "DropletUtils", "Matrix", "Seurat", "SeuratObject",
  "SingleCellExperiment", "SingleR", "SoupX", "celldex", "edgeR", "fgsea",
  "ggplot2", "ggrepel", "gypsum", "harmony", "jsonlite", "liana", "limma",
  "lme4", "nnls", "patchwork", "readxl", "scDblFinder", "scales", "scater",
  "scuttle", "svglite", "ragg", "yaml", "dplyr", "future", "statmod", "metafor"
)
available <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
if (!all(available)) {
  stop("Missing required R packages: ", paste(required[!available], collapse = ", "),
       ". Create environment.yml and run scripts/00_install_locked_r_dependencies.R")
}
if (!identical(utils::packageDescription("liana")$RemoteSha,
               "6cab46c54234f861ea176c3de77c4b8aa45ecb3d")) {
  stop("LIANA source commit mismatch; run scripts/00_install_locked_r_dependencies.R")
}
cat("Environment preflight PASS:", length(required), "required R packages available.\n")
if ("--original-runtime" %in% commandArgs(TRUE)) {
  original <- read.delim("environment.lock.tsv", stringsAsFactors = FALSE)
  observed <- vapply(original$component, function(package) {
    if (package == "R") return(as.character(getRversion()))
    if (!requireNamespace(package, quietly = TRUE)) return("MISSING")
    as.character(packageVersion(package))
  }, character(1))
  mismatch <- observed != original$version
  if (any(mismatch)) stop(
    "Original-data reproduction requires the recorded runtime versions: ",
    paste(paste0(original$component[mismatch], " expected ", original$version[mismatch],
                 ", found ", observed[mismatch]), collapse = "; "),
    ". The supplied Linux bootstrap is certified for retained-result checks; ",
    "restore environment.lock.tsv versions before original discovery processing.")
  cat("Original-runtime version preflight PASS.\n")
}
