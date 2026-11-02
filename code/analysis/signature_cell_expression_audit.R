#!/usr/bin/env Rscript
# Run from repository root. Retained singlets, same lineage labels in both assays.
# No cell-level hypothesis test is performed.
suppressPackageStartupMessages({library(Seurat);library(Matrix)})
args<-commandArgs(trailingOnly=TRUE)
input<-if(length(args))args[1] else "data/processed/scrna/scrna_annotated.rds"
source("scripts/lib_common.R")
out<-if(length(args)>=2) args[2] else "analysis"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_tsv<-function(x,n) write.table(x,file.path(out,paste0("signature_",n,".tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
obj<-readRDS(input)
meta<-obj@meta.data
stopifnot(nrow(meta)==14846, all(c("sample_id","patient_id","side","cell_type") %in% names(meta)),
          all(c("RNA","SoupX") %in% Assays(obj)), !anyNA(meta[,c("sample_id","patient_id","side","cell_type")]),
          length(unique(meta$sample_id))==4, length(unique(meta$cell_type))==9)
validate_paired_samples(unique(meta[,c("sample_id","patient_id","side")]),expected_pairs=2)
sig<-read.delim("results/tables/bulk/frozen_multicellular_signature.tsv",check.names=FALSE)
target<-unique(c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5",
 "CKM","MYH1","MYH2","MYH7","TNNT3","TTN","H19","ADIPOQ","TENT5A","ESR1","METTL3","MSTN","SMAD3",
 sig$gene))
availability<-do.call(rbind,lapply(c("RNA","SoupX"),function(a) data.frame(assay=a,gene=target,
 available=target %in% rownames(obj[[a]]))))
write_tsv(availability,"requested_gene_availability")
# Absent genes are explicit in the manifest and are not imputed as zero counts.
target<-Reduce(intersect,list(target,rownames(obj[["RNA"]]),rownames(obj[["SoupX"]])))
stopifnot(length(target)>0)
libraries<-sort(unique(meta$sample_id));lineages<-sort(unique(meta$cell_type))
audit<-list();cell_output<-list();assay_manifest<-list()
for(assay in c("RNA","SoupX")) {
 message("Audit ",assay)
 # Join retained library-specific layers without changing counts or labels.
 obj<-JoinLayers(obj,assay=assay)
 m<-LayerData(obj,assay=assay,layer="counts")
 m<-m[,rownames(meta),drop=FALSE]
 totals<-Matrix::colSums(m)
 selected<-as.matrix(m[target,,drop=FALSE])
 cp10k<-sweep(selected,2,pmax(totals,1),"/")*1e4
 log_norm<-log1p(cp10k)
 assay_manifest[[assay]]<-data.frame(assay=assay,stage=if(assay=="RNA")"before_ambient_correction" else "after_SoupX",
 n_genes=nrow(m),n_cells=ncol(m),n_target_genes=nrow(selected),
 transformation="UMI_counts;per_cell_counts_divided_by_all_gene_UMIs_times_10000;natural_log1p",source_object=input)
 for(lib in libraries)for(ct in lineages) {
  ix<-which(meta$sample_id==lib&meta$cell_type==ct)
  if(!length(ix)) {
   audit[[paste(assay,lib,ct)]]<-data.frame(assay=assay,sample_id=lib,patient_id=sub("__.*$","",lib),side=sub("^.*__","",lib),
    cell_type=ct,gene=target,n_cells=0,n_detected=0,detection_fraction=NA,mean_UMI_per_cell=NA,mean_UMI_detected_cells=NA,
    mean_cp10k=NA,mean_log1p_cp10k=NA,median_log1p_cp10k=NA,p90_log1p_cp10k=NA,max_log1p_cp10k=NA,
    inferential_status="descriptive_library_lineage_distribution_no_cell_independence")
   next
  }
  counts<-selected[,ix,drop=FALSE];ln<-log_norm[,ix,drop=FALSE];det<-rowSums(counts>0)
  audit[[paste(assay,lib,ct)]]<-data.frame(assay=assay,sample_id=lib,patient_id=meta$patient_id[ix[1]],side=meta$side[ix[1]],
   cell_type=ct,gene=target,n_cells=length(ix),n_detected=det,detection_fraction=det/length(ix),mean_UMI_per_cell=rowMeans(counts),
   mean_UMI_detected_cells=ifelse(det>0,rowSums(counts)/det,NA),mean_cp10k=rowMeans(cp10k[,ix,drop=FALSE]),
   mean_log1p_cp10k=rowMeans(ln),median_log1p_cp10k=apply(ln,1,median),p90_log1p_cp10k=apply(ln,1,quantile,.9),
   max_log1p_cp10k=apply(ln,1,max),inferential_status="descriptive_library_lineage_distribution_no_cell_independence")
 }
 # Provide per-cell observations for unexpected/nominated genes as source data.
 unexpected<-intersect(c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5","ESR1"),target)
 cell_output[[assay]]<-do.call(rbind,lapply(unexpected,function(g) data.frame(assay=assay,cell_id=rownames(meta),sample_id=meta$sample_id,
  patient_id=meta$patient_id,side=meta$side,cell_type=meta$cell_type,gene=g,UMI=selected[g,],log1p_cp10k=log_norm[g,])))
 rm(m,selected,cp10k,log_norm);gc(verbose=FALSE)
}
audit<-do.call(rbind,audit);write_tsv(audit,"expression_by_library_lineage_before_after")
write_tsv(do.call(rbind,assay_manifest),"expression_assay_manifest")
connection<-gzfile(file.path(out,"signature_unexpected_gene_cell_source.tsv.gz"),"wt")
write.table(do.call(rbind,cell_output),connection,sep="\t",quote=FALSE,row.names=FALSE,na="NA");close(connection)
# Join the two correction stages for direct descriptive changes.
before<-audit[audit$assay=="RNA",];after<-audit[audit$assay=="SoupX",]
key<-c("sample_id","patient_id","side","cell_type","gene","n_cells")
comparison<-merge(before[,c(key,"detection_fraction","mean_UMI_per_cell","mean_log1p_cp10k")],
 after[,c(key,"detection_fraction","mean_UMI_per_cell","mean_log1p_cp10k")],by=key,suffixes=c(".before",".after"))
comparison$detection_fraction_change<-comparison$detection_fraction.after-comparison$detection_fraction.before
comparison$mean_UMI_change<-comparison$mean_UMI_per_cell.after-comparison$mean_UMI_per_cell.before
write_tsv(comparison,"expression_correction_changes")
# Expression magnitudes for every displayed signature gene and its chosen lineage.
displayed<-merge(sig,audit[audit$assay=="SoupX",],by=c("gene","cell_type"),all.x=TRUE)
write_tsv(displayed,"frozen_gene_expression_all_four_libraries")
run_manifest<-data.frame(
 field=c("input_object","input_object_md5","audit_script_md5","signature_md5",
         "R_version","Seurat_version","n_retained_singlets","n_libraries","n_lineages",
         "n_requested_genes","n_available_genes_in_both_assays","n_summary_rows",
         "scientific_scope"),
 value=c(input,unname(tools::md5sum(input)),
         unname(tools::md5sum("analysis/signature_cell_expression_audit.R")),
         unname(tools::md5sum("results/tables/bulk/frozen_multicellular_signature.tsv")),
         as.character(getRversion()),as.character(packageVersion("Seurat")),
         nrow(meta),length(libraries),length(lineages),nrow(availability)/2,
         length(target),nrow(audit),"descriptive_per_cell_before_after_same_retained_singlets_no_cell_level_hypothesis_test"))
write_tsv(run_manifest,"cell_expression_run_manifest")
message("Cell-expression audit complete: ",nrow(audit)," gene x lineage x library x assay summaries.")
