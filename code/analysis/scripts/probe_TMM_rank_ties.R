#!/usr/bin/env Rscript
# Numerical amendment probe only; does not run biological tests or change primary factors.
# Adapted from edgeR 4.4.2 .calcFactorTMM/calcNormFactors.default.
suppressPackageStartupMessages(library(edgeR))
argv <- commandArgs(trailingOnly = TRUE)
allowed <- c("--input-dir", "--output-dir", "--platform")
keys <- sub("=.*$", "", argv)
if (any(!keys %in% allowed) || anyDuplicated(keys) || any(!grepl("=.+$", argv))) stop("Unknown/duplicate/empty argument")
arg <- function(k, fallback) {x <- argv[startsWith(argv, paste0(k, "="))]; if(length(x)) sub(paste0("^",k,"="), "", x) else fallback}
inputs <- arg("--input-dir", "analysis/bulk_gate_inputs")
out <- arg("--output-dir", "logs/TMM_rank_probe")
platform <- arg("--platform", Sys.info()[["sysname"]])
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(digits = 17)
rd <- function(n) read.delim(if(grepl("gz$",n)) gzfile(file.path(inputs,n)) else file.path(inputs,n),check.names=FALSE)
wr <- function(x,n) write.table(x,file.path(out,n),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
tab <- rd("GSE311180_tested_counts.tsv.gz"); counts <- as.matrix(tab[,-1]);rownames(counts)<-tab$gene
meta <- rd("samples.tsv"); meta <- meta[meta$cohort=="GSE311180",]
gates <- rd("RNA_expected_pairs.tsv"); keep<-sort(gates$patient_id[gates$interpretation=="consistent_RNA_profiles"])
stopifnot(identical(keep,c("STB-050","STB-062","STB-085","STB-086","STB-095")),
          identical(colnames(counts),meta$sample_id),nrow(counts)==17961L)
counts <- counts[, meta$patient_id %in% keep, drop=FALSE]
libsize <- colSums(counts)
x <- counts[rowSums(counts>0)>0,,drop=FALSE]
f75 <- suppressWarnings(getFromNamespace(".calcFactorQuantile","edgeR")(data=x,lib.size=libsize,p=.75))
refColumn <- if(median(f75)<1e-20) which.max(colSums(sqrt(x))) else which.min(abs(f75-mean(f75)))
rank_digits <- 12L; factor_digits <- 14L
changed <- list(); details <- list()
one <- function(obs,ref,nO,nR,sample_id) {
  logR <- log2((obs/nO)/(ref/nR)); absE <- (log2(obs/nO)+log2(ref/nR))/2
  v <- (nO-obs)/nO/obs+(nR-ref)/nR/ref
  fin <- is.finite(logR)&is.finite(absE)&absE> -1e10
  genes <- rownames(x)[fin];logR<-logR[fin];absE<-absE[fin];v<-v[fin]
  # The same reference-library early return as native edgeR.
  if(max(abs(logR))<1e-6) {
    details[[sample_id]] <<- data.frame(sample_id=sample_id,n_finite=length(logR),changed_M_ranks=0,
      changed_A_ranks=0,changed_retained_genes=0,native_retained=NA,stabilized_retained=NA)
    return(c(native=1,stabilized=1))
  }
  n<-length(logR);loL<-floor(n*.3)+1;hiL<-n+1-loL;loS<-floor(n*.05)+1;hiS<-n+1-loS
  rM<-rank(logR);rA<-rank(absE)
  stableM<-rank(round(logR,rank_digits));stableA<-rank(round(absE,rank_digits))
  native_keep<-(rM>=loL&rM<=hiL)&(rA>=loS&rA<=hiS)
  stable_keep<-(stableM>=loL&stableM<=hiL)&(stableA>=loS&stableA<=hiS)
  weighted<-function(z) {f<-sum(logR[z]/v[z],na.rm=TRUE)/sum(1/v[z],na.rm=TRUE);if(is.na(f))f<-0;2^f}
  native<-weighted(native_keep);stable<-weighted(stable_keep)
  expected<-getFromNamespace(".calcFactorTMM","edgeR")(obs,ref,libsize.obs=nO,libsize.ref=nR,
    logratioTrim=.3,sumTrim=.05,doWeighting=TRUE,Acutoff=-1e10)
  stopifnot(abs(native-expected)<1e-14)
  details[[sample_id]] <<- data.frame(sample_id=sample_id,n_finite=n,changed_M_ranks=sum(rM!=stableM),
    changed_A_ranks=sum(rA!=stableA),changed_retained_genes=sum(native_keep!=stable_keep),
    native_retained=sum(native_keep),stabilized_retained=sum(stable_keep))
  if(any(native_keep!=stable_keep)) changed[[sample_id]] <<- data.frame(sample_id=sample_id,
    gene=genes[native_keep!=stable_keep],logR=logR[native_keep!=stable_keep],absE=absE[native_keep!=stable_keep],
    native_keep=native_keep[native_keep!=stable_keep],stable_keep=stable_keep[native_keep!=stable_keep])
  c(native=native,stabilized=stable)
}
raw <- vapply(seq_len(ncol(x)),function(j) one(x[,j],x[,refColumn],libsize[j],libsize[refColumn],colnames(x)[j]),numeric(2))
native<-raw[1,]/exp(mean(log(raw[1,])));stable<-raw[2,]/exp(mean(log(raw[2,])))
canonical<-round(stable,factor_digits)
builtin<-calcNormFactors(counts,method="TMM")
stopifnot(max(abs(native-builtin))<1e-14,max(abs(canonical-stable))<=5.1e-15)
factor_table<-data.frame(sample_id=colnames(counts),reference_sample=colnames(counts)[refColumn],
  native_TMM_factor=as.numeric(builtin),stabilized_rank_TMM_factor=as.numeric(stable),
  canonical_stabilized_rank_TMM_factor=as.numeric(canonical),native_minus_stabilized=native-stable,
  rank_decimal_places=rank_digits,factor_decimal_places=factor_digits,platform=platform)
wr(factor_table,"normalization_probe.tsv")
wr(data.frame(sample_id=colnames(counts),native_hex=sprintf("%a",builtin),
  stabilized_hex=sprintf("%a",stable),canonical_hex=sprintf("%a",canonical)),"factor_double_hex.tsv")
wr(do.call(rbind,details),"trim_mask_summary.tsv")
if(length(changed))wr(do.call(rbind,changed),"changed_trim_genes.tsv") else wr(data.frame(
  sample_id=character(),gene=character(),logR=numeric(),absE=numeric(),native_keep=logical(),stable_keep=logical()),"changed_trim_genes.tsv")
wr(data.frame(sample_id=colnames(counts),factor=format(canonical,digits=17,scientific=FALSE,trim=TRUE)),"canonical_factors.tsv")
wr(data.frame(package=c("R","edgeR"),version=c(as.character(getRversion()),as.character(packageVersion("edgeR")))),"runtime_versions.tsv")
manifest<-list(platform=platform,rank_decimal_places=rank_digits,factor_decimal_places=factor_digits,
  counts_sha256=digest::digest(file=file.path(inputs,"GSE311180_tested_counts.tsv.gz"),algo="sha256"),
  script_sha256=digest::digest(file=sub("^--file=","",commandArgs()[startsWith(commandArgs(),"--file=")][1]),algo="sha256"),
  native_TMM_source=digest::digest(paste(deparse(getFromNamespace(".calcFactorTMM","edgeR")),collapse="\n"),algo="sha256"),
  native_copy_max_absolute_error=max(abs(native-builtin)),
  max_native_vs_stabilized=max(abs(native-stable)),max_canonical_rounding=max(abs(canonical-stable)),
  primary_archived_factors_modified=FALSE,biological_tests_run=FALSE)
jsonlite::write_json(manifest,file.path(out,"manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=17)
cat("NUMERICAL_TMM_PROBE_COMPLETED",platform,"native/stable",max(abs(native-stable)),"\n")
