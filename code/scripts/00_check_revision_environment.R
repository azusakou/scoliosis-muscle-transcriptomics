#!/usr/bin/env Rscript
packages <- c("edgeR", "statmod", "metafor", "lme4", "ggplot2", "patchwork", "ggrepel", "scales", "svglite", "ragg")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing revision dependencies: ", paste(missing, collapse = ", "),
                          ". Install the dependencies in environment.yml before reproduction.")
if (as.character(packageVersion("metafor")) != "4.8.0")
  stop("Retained revision reproduction requires metafor 4.8.0; found ", as.character(packageVersion("metafor")))
cat("Revision dependency preflight PASS.\n")
