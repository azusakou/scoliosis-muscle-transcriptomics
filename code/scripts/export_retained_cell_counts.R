#!/usr/bin/env Rscript
# Export genuine fixed-target sparse counts only after exact retained-source joins.
# Usage: Rscript scripts/export_retained_cell_counts.R object.rds [output_directory]
suppressPackageStartupMessages({library(Seurat);library(Matrix)})
argv <- commandArgs(trailingOnly=TRUE)
if(!length(argv))stop("Supply a reconstructed annotated Seurat RDS; no object is bundled")
script <- normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]),mustWork=TRUE)
root <- dirname(dirname(script))
input <- normalizePath(argv[1],mustWork=TRUE)
out <- if(length(argv)>=2)argv[2] else file.path(root,"data/retained")
hash <- function(p)digest::digest(file=p,algo="sha256")
read_tsv <- function(p)read.delim(if(grepl("\\.gz$",p))gzfile(p) else p,check.names=FALSE,stringsAsFactors=FALSE)
write_tsv <- function(x,p) {
 con <- if(grepl("\\.gz$",p))gzfile(p,"wt") else file(p,"wt")
 on.exit(close(con));write.table(x,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
required <- c("sample_id","patient_id","side","cell_type")
metadata_file <- file.path(root,"results/tables/annotation/annotated_cell_metadata.tsv.gz")
metadata <- read_tsv(metadata_file)
stopifnot(nrow(metadata)==14846,!anyDuplicated(metadata$cell_id))
sig <- read_tsv(file.path(root,"results/tables/bulk/frozen_multicellular_signature.tsv"))
requested <- unique(c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5",
 "CKM","MYH1","MYH2","MYH7","TNNT3","TTN","H19","ADIPOQ","TENT5A","ESR1","METTL3","MSTN","SMAD3",sig$gene))
biology <- c("ANGPT1","LEP","ADIPOQ","PDGFRA","COL1A1","COL1A2","COL3A1","DCN","LUM","PI16","DPP4",
 "CFD","APOD","POSTN","IGKC","ACTA1","MYL2","DES")
targets <- unique(c(requested,biology))
stopifnot(length(requested)==158,length(biology)==18)
object <- readRDS(input); md <- object@meta.data
stopifnot(nrow(md)==nrow(metadata),!anyDuplicated(rownames(md)),setequal(rownames(md),metadata$cell_id))
md <- md[match(metadata$cell_id,rownames(md)),,drop=FALSE]
labels <- c(required,"seurat_clusters","cell_state_label","annotation_confidence","passes_qc_metrics","passes_qc","scDblFinder.class")
for(f in labels)if(!identical(as.character(md[[f]]),as.character(metadata[[f]])))stop("Object/retained metadata label mismatch: ",f)
stopifnot(length(unique(md$sample_id))==4,length(unique(md$patient_id))==2,length(unique(md$cell_type))==9)
# Stored decimal QC fields retain their existing serialization-level criterion.
qc <- c("subsets_Mito_percent","subsets_Ribo_percent","scDblFinder.score")
for(f in qc)if(any(abs(as.numeric(md[[f]])-metadata[[f]])>2e-12*pmax(1,abs(metadata[[f]]))))stop("Object/retained QC mismatch: ",f)
sources <- c(unexpected="signature_unexpected_gene_cell_source.tsv.gz",biology="signature_biology_all_target_cell_source.tsv.gz")
source_tables <- lapply(sources,function(n)read_tsv(file.path(root,"analysis",n)))
expected_source_rows <- c(unexpected=296920L,biology=534456L)
expected_source_genes <- list(unexpected=c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5","ESR1"),biology=biology)
for(n in names(sources)) {
 z<-source_tables[[n]]
 stopifnot(nrow(z)==expected_source_rows[[n]],!anyDuplicated(z[,c("assay","cell_id","gene")]),
  setequal(z$gene,expected_source_genes[[n]]),setequal(z$assay,c("RNA","SoupX")),all(z$cell_id%in%metadata$cell_id))
 for(f in required)if(!identical(as.character(z[[f]]),as.character(metadata[[f]][match(z$cell_id,metadata$cell_id)])))stop("Existing source label mismatch: ",n,"/",f)
 stopifnot(all(table(z$assay,z$gene)==nrow(metadata)))
}
counts <- list();features <- list();availability <- list();assays <- list()
for(a in c("RNA","SoupX")) {
 message("Source match ",a)
 object <- JoinLayers(object,assay=a)
 m <- LayerData(object,assay=a,layer="counts")[,metadata$cell_id,drop=FALSE]
 stopifnot(!anyDuplicated(rownames(m)),ncol(m)==nrow(metadata),
  all(is.finite(m@x)),all(m@x>=0),all(m@x==round(m@x)))
 totals <- Matrix::colSums(m); detected <- Matrix::colSums(m>0)
 for(prefix in c("nCount_","nFeature_")) {
  actual<-if(prefix=="nCount_")totals else detected;field<-paste0(prefix,a)
  if(!all(actual==metadata[[field]])||!all(actual==md[[field]]))stop("Exact full-assay count/feature total mismatch: ",field)
 }
 selected <- m[intersect(targets,rownames(m)),,drop=FALSE]
 for(n in names(sources)) {
  z<-source_tables[[n]];z<-z[z$assay==a,]
  stopifnot(all(z$gene%in%rownames(selected)))
  actual<-as.numeric(selected[cbind(match(z$gene,rownames(selected)),match(z$cell_id,metadata$cell_id))])
  if(!all(actual==z$UMI))stop("Exact existing UMI mismatch: ",n,"/",a)
  fulltotals<-as.numeric(totals[match(z$cell_id,metadata$cell_id)])
  if("all_gene_total_UMI"%in%names(z)&&!all(fulltotals==z$all_gene_total_UMI))stop("Existing all-gene total mismatch")
  normalized<-actual/fulltotals*1e4
  if(any(abs(log1p(normalized)-z$log1p_cp10k)>2e-10*pmax(1,abs(z$log1p_cp10k))))stop("Existing normalization mismatch")
 }
 nonzero<-summary(selected);nonzero<-nonzero[nonzero$x>0,]
 counts[[a]]<-data.frame(assay=a,cell_id=colnames(selected)[nonzero$j],gene=rownames(selected)[nonzero$i],UMI=nonzero$x)
 features[[a]]<-data.frame(assay=a,gene=rownames(m))
 availability[[a]]<-data.frame(assay=a,gene=targets,available=targets%in%rownames(m),
  original_expression_requested=targets%in%requested,closed_biology_target=targets%in%biology)
 assays[[a]]<-data.frame(assay=a,n_genes=nrow(m),n_cells=ncol(m),n_export_target_genes=nrow(selected),
  n_export_nonzero_records=nrow(nonzero),all_gene_total_source_field=paste0("nCount_",a))
 rm(m,selected);gc(verbose=FALSE)
}
# Validate the earlier availability record too; an unavailable gene is never zero-filled.
old_availability<-read_tsv(file.path(root,"analysis/signature_requested_gene_availability.tsv"))
new_availability<-do.call(rbind,availability)
for(i in seq_len(nrow(old_availability))) {
 z<-old_availability[i,];actual<-new_availability$available[new_availability$assay==z$assay&new_availability$gene==z$gene]
 stopifnot(length(actual)==1,identical(actual,z$available))
}
names_out<-c("cell_target_counts.tsv.gz","cell_assay_features.tsv.gz","cell_target_availability.tsv","cell_assay_manifest.tsv","retained_cell_counts_export.json")
if(any(file.exists(file.path(out,names_out))))stop("Output exists; use a separate new export directory")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_tsv(do.call(rbind,counts),file.path(out,names_out[1]))
write_tsv(do.call(rbind,features),file.path(out,names_out[2]))
write_tsv(new_availability,file.path(out,names_out[3]))
write_tsv(do.call(rbind,assays),file.path(out,names_out[4]))
files<-lapply(names_out[1:4],function(n)list(path=n,bytes=unname(file.info(file.path(out,n))$size),sha256=hash(file.path(out,n))))
report<-list(status="PASS",scope="Actual reconstructed-object count export; no cell-level inference, no source relabelling or imputed values",
 input_object_role="independently_reconstructed_public_annotated_singlets",input_object_sha256=hash(input),
 exporter_sha256=hash(script),metadata_relative_path="results/tables/annotation/annotated_cell_metadata.tsv.gz",
 metadata_sha256=hash(metadata_file),signature_sha256=hash(file.path(root,"results/tables/bulk/frozen_multicellular_signature.tsv")),
 source_rows_exactly_matched=as.list(expected_source_rows),source_sha256=lapply(sources,function(n)hash(file.path(root,"analysis",n))),
 n_cells=nrow(metadata),n_expression_requested_genes=length(requested),n_biology_targets=length(biology),
 n_union_targets=length(targets),assay_manifest=do.call(rbind,assays),files=files,
 zero_encoding="Only available assay/gene combinations may interpret an omitted cell count as zero; unavailable genes remain missing",
 total_UMI_scope="Complete assay UMI sums exactly checked against retained metadata; not target-subset sums",
 R_version=as.character(getRversion()),Seurat_version=as.character(packageVersion("Seurat")))
jsonlite::write_json(report,file.path(out,names_out[5]),pretty=TRUE,auto_unbox=TRUE,digits=NA)
message("Export PASS: ",length(targets)," fixed union targets; ",sum(vapply(counts,nrow,integer(1)))," nonzero records")
