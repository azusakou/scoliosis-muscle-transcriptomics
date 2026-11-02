#!/usr/bin/env Rscript
# Minimal nonduplicated counts for the two original script61 gene-level models.
suppressPackageStartupMessages({library(readxl);library(edgeR)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1L || length(args)>2L) stop("Usage: Rscript scripts/export_specificity_counts.R GSE254300_readcount.xlsx [NEW_OUTPUT_DIRECTORY]")
script <- normalizePath(sub("^--file=","",commandArgs(FALSE)[grep("^--file=",commandArgs(FALSE))][1]),mustWork=TRUE)
root <- dirname(dirname(script)); input <- normalizePath(args[1],mustWork=TRUE)
out <- if(length(args)>1L)args[2] else file.path(root,"data/retained/gse254300_specificity")
if(file.exists(out)) stop("Specificity count destination already exists; use a new directory")
hash <- function(p)digest::digest(file=p,algo="sha256")
rd <- function(p)read.delim(if(grepl("[.]gz$",p))gzfile(p) else p,check.names=FALSE,stringsAsFactors=FALSE)
matrix_read <- function(p){x<-rd(p);m<-as.matrix(x[,-1,drop=FALSE]);rownames(m)<-x[[1]];storage.mode(m)<-"double";m}
metadata_path <- file.path(root,"metadata/audit/bulk_samples.tsv")
cs_path <- file.path(root,"results/tables/disease_control_gse254300_cs/GSE254300_CS_paired_edgeR_full.tsv.gz")
joint_path <- file.path(root,"results/tables/disease_control_gse254300_cs/GSE254300_AIS_by_side_interaction_edgeR_full.tsv.gz")
a_path <- file.path(root,"analysis/bulk_gate_inputs/GSE254300_tested_counts.tsv.gz")
audit <- rd(metadata_path); audit <- audit[audit$cohort %in% c("AIS","congenital_scoliosis"),]
stopifnot(nrow(audit)==20L,all(table(audit$cohort)==10L),all(table(audit$cohort,audit$side)==5L),all(table(audit$patient_id)==2L))
raw <- as.data.frame(read_excel(input));raw$GeneID <- make.unique(as.character(raw$GeneID))
stopifnot(!anyDuplicated(raw$GeneID),all(audit$processed_sample_id %in% names(raw)))
counts <- as.matrix(raw[,audit$processed_sample_id,drop=FALSE]);rownames(counts)<-raw$GeneID;storage.mode(counts)<-"double"
stopifnot(all(is.finite(counts)),all(counts>=0),all(counts==round(counts)))
info <- audit;info$cohort_label<-ifelse(info$cohort=="AIS","AIS","CS");info$side<-factor(info$side,levels=c("convex","concave"))
cs <- info$cohort_label=="CS"
keep_joint <- filterByExpr(counts,group=interaction(info$cohort_label,info$side),min.count=10)
keep_cs <- filterByExpr(counts[,cs,drop=FALSE],group=info$side[cs],min.count=10)
genes_joint <- rd(joint_path)$gene;genes_cs<-rd(cs_path)$gene
stopifnot(length(genes_joint)==14967L,length(genes_cs)==14450L,
 identical(rownames(counts)[keep_joint],genes_joint),identical(rownames(counts)[keep_cs],genes_cs),
 all(genes_cs %in% genes_joint))
A <- matrix_read(a_path);stopifnot(nrow(A)==14296L,ncol(A)==10L,all(rownames(A)%in%genes_joint),
 identical(colnames(A),info$processed_sample_id[!cs]),
 identical(unname(counts[rownames(A),colnames(A),drop=FALSE]),unname(A)))
extra <- genes_joint[!genes_joint%in%rownames(A)];stopifnot(length(extra)==671L)
write_matrix <- function(m,p){con<-gzfile(p,"wt");on.exit(close(con));write.table(data.frame(gene=rownames(m),m,check.names=FALSE),con,sep="\t",quote=FALSE,row.names=FALSE)}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
names <- c("CS_joint_universe_counts.tsv.gz","AIS_joint_extra_counts.tsv.gz")
write_matrix(counts[genes_joint,cs,drop=FALSE],file.path(out,names[1]))
write_matrix(counts[extra,!cs,drop=FALSE],file.path(out,names[2]))
report <- list(status="PASS",scope="Exact original public XLSX counts and original script61 tested universes; existing 14296-gene AIS counts reused without duplication",
 source_public_role="GSE254300 deposited count XLSX",source_basename=basename(input),source_sha256=hash(input),
 exporter_sha256=hash(script),original_model_producer_sha256=hash(file.path(root,"scripts/61_gse254300_cs_specificity.R")),
 metadata_sha256=hash(metadata_path),existing_A_counts_sha256=hash(a_path),
 original_CS_result_sha256=hash(cs_path),original_interaction_result_sha256=hash(joint_path),
 original_full_matrix_genes=nrow(counts),n_CS_tested_genes=length(genes_cs),n_joint_tested_genes=length(genes_joint),
 n_AIS_reused_genes=nrow(A),n_AIS_extra_genes=length(extra),n_AIS_libraries=sum(!cs),n_CS_libraries=sum(cs),
 full_matrix_filterByExpr_verified=TRUE,raw_AIS_reused_count_values_exactly_equal=TRUE,
 filtering=list(joint="filterByExpr(counts_all, group=interaction(cohort_label,side), min.count=10)",CS="filterByExpr(counts_cs, group=side, min.count=10)"),
 files=lapply(names,function(n)list(path=n,bytes=unname(file.info(file.path(out,n))$size),sha256=hash(file.path(out,n)))),
 R_version=as.character(getRversion()),edgeR_version=as.character(packageVersion("edgeR")))
jsonlite::write_json(report,file.path(out,"specificity_counts_export.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
message("PASS: exported CS14967x10 + AIS671x10; original AIS14296x10 reused and both original full-matrix gene filters exactly verified")
