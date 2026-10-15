#!/usr/bin/env Rscript

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("Install r-remotes from environment.yml before running this script")
}
# Conda/Bioconda can resolve older releases within Bioconductor 3.20. Install
# the original statistical versions explicitly rather than silently changing
# QL estimation when the environment is rebuilt.
locked <- data.frame(
  package = c("GenomeInfoDbData", "celldex", "statmod", "limma", "edgeR", "basilisk.utils", "basilisk"),
  version = c("1.2.13", "1.16.0", "1.5.1", "3.62.2", "4.4.2", "1.18.0", "1.18.0"),
  url = c("https://bioconductor.org/packages/3.20/data/annotation/src/contrib/GenomeInfoDbData_1.2.13.tar.gz",
          "https://bioconductor.org/packages/3.20/data/experiment/src/contrib/celldex_1.16.0.tar.gz",
          "https://cran.r-project.org/src/contrib/Archive/statmod/statmod_1.5.1.tar.gz",
          "https://bioconductor.org/packages/3.20/bioc/src/contrib/limma_3.62.2.tar.gz",
          "https://bioconductor.org/packages/3.20/bioc/src/contrib/edgeR_4.4.2.tar.gz",
          "https://bioconductor.org/packages/3.20/bioc/src/contrib/basilisk.utils_1.18.0.tar.gz",
          "https://bioconductor.org/packages/3.20/bioc/src/contrib/basilisk_1.18.0.tar.gz"),
  md5 = c("8798e58b89d84d9fa6a85f2191704668", "ba3dc1f475b673ee0303293d6d997194",
          "66477fbc8734e31c5985aff1a57b9a7e", "76032128ad85285ed3dfc05a73cb9b51",
          "4ce1d994f76cc821f90d86e09fef6096", "89953864457a1974b657ed7dc01f8d50",
          "d75135397396224a3d7e2fe691aadb33"))
for (i in seq_len(nrow(locked))) {
  package <- locked$package[i]
  matches <- requireNamespace(package, quietly = TRUE) &&
    identical(as.character(packageVersion(package)), locked$version[i])
  if (!matches) {
    archive <- tempfile(fileext = ".tar.gz")
    # Bioconda's data-package post-link hooks are order/PATH dependent in an
    # explicit-lock installation. Install these sources explicitly as well.
    # The archive mirror is accepted only with the same recorded MD5.
    urls <- locked$url[i]
    if (package %in% c("GenomeInfoDbData", "celldex")) {
      urls <- c(urls, paste0("https://bioarchive.galaxyproject.org/", basename(urls)))
    }
    verified <- FALSE
    for (url in urls) {
      downloaded <- tryCatch({
        download.file(url, archive, mode = "wb", method = "libcurl")
        file.exists(archive) && unname(tools::md5sum(archive)) == locked$md5[i]
      }, error = function(e) FALSE)
      if (isTRUE(downloaded)) { verified <- TRUE; break }
    }
    if (!verified) stop("No checksum-verified source archive for ", package)
    # Do not load a prior version before installing a replacement in this R session.
    if (package %in% loadedNamespaces()) unloadNamespace(package)
    install.packages(archive, repos = NULL, type = "source")
    unlink(archive)
  }
}
# LIANA declares floating GitHub Remotes for basilisk, even though this analysis
# uses Bioconductor3.20. Install its hard dependencies via environment.yml, then
# suppress automatic Remotes processing when installing the fixed LIANA commit.
omnipath_commit <- "3e1139172d3036ffd772ebb945bc5e338ac63fde"
omnipath_matches <- requireNamespace("OmnipathR", quietly = TRUE) &&
  identical(as.character(packageVersion("OmnipathR")), "4.1.0") &&
  identical(packageDescription("OmnipathR")$RemoteSha, omnipath_commit)
if (!omnipath_matches) {
  if ("OmnipathR" %in% loadedNamespaces()) unloadNamespace("OmnipathR")
  remotes::install_github(paste0("saezlab/OmnipathR@", omnipath_commit),
                          upgrade = "never", dependencies = NA)
}
stopifnot(as.character(packageVersion("OmnipathR")) == "4.1.0",
          packageDescription("OmnipathR")$RemoteSha == omnipath_commit)

target_version <- "0.1.14"
target_commit <- "6cab46c54234f861ea176c3de77c4b8aa45ecb3d"
installed_matches <- function() {
  if (!requireNamespace("liana", quietly = TRUE)) return(FALSE)
  description <- utils::packageDescription("liana")
  identical(as.character(utils::packageVersion("liana")), target_version) &&
    identical(description$RemoteSha, target_commit)
}
needs_liana <- !installed_matches()
if (needs_liana) {
  remotes::install_github(paste0("saezlab/liana@", target_commit), upgrade = "never",
                          dependencies = FALSE)
}
# Source overlays can replace dependencies already loaded in this installer
# session. Check the final namespace graph in a fresh R process, which is also
# how the scientific workflow starts, rather than accepting a stale namespace.
verification <- paste0(
  'stopifnot(requireNamespace("liana", quietly=TRUE), ',
  'as.character(packageVersion("liana")) == "0.1.14", ',
  'packageDescription("liana")$RemoteSha == "', target_commit, '", ',
  'requireNamespace("celldex", quietly=TRUE), ',
  'requireNamespace("GenomeInfoDbData", quietly=TRUE)); ',
  'cat("Fresh-process locked namespace verification PASS\\n")')
status <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(verification)))
stopifnot(status == 0L)
for (i in seq_len(nrow(locked))) {
  stopifnot(as.character(packageVersion(locked$package[i])) == locked$version[i])
}
cat("Locked non-conda dependency PASS: liana", target_version,
    "at source commit 6cab46c54234f861ea176c3de77c4b8aa45ecb3d.\n")
