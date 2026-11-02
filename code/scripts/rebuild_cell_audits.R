#!/usr/bin/env Rscript
# Recompute the full fixed-target audits from genuine retained sparse counts.
# Usage: Rscript scripts/rebuild_cell_audits.R [counts_directory] [new_run_directory]
# No Seurat, full annotated object, old summary input, cell test or relabelling.
suppressPackageStartupMessages(library(Matrix))
args<-commandArgs(trailingOnly=TRUE)
script<-normalizePath(sub("^--file=","",commandArgs(FALSE)[grep("^--file=",commandArgs(FALSE))][1]),mustWork=TRUE)
root<-dirname(dirname(script))
input<-if(length(args))normalizePath(args[1],mustWork=TRUE) else file.path(root,"data/retained")
logs<-file.path(root,"logs/cell_audit_replay");dir.create(logs,recursive=TRUE,showWarnings=FALSE)
run<-if(length(args)>=2)args[2] else tempfile("run-",tmpdir=logs)
if(file.exists(run))stop("Replay output already exists; use a separate new run directory")
hash<-function(p)digest::digest(file=p,algo="sha256")
read_tsv<-function(p)read.delim(if(grepl("\\.gz$",p))gzfile(p) else p,check.names=FALSE,stringsAsFactors=FALSE)
metadata_path<-file.path(root,"results/tables/annotation/annotated_cell_metadata.tsv.gz")
binding_path<-if(length(args)&&file.exists(file.path(input,"retained_cell_counts_export.json")))
 file.path(input,"retained_cell_counts_export.json") else file.path(root,"metadata/audit/retained_cell_counts_export.json")
binding<-jsonlite::fromJSON(binding_path,simplifyVector=FALSE)
stopifnot(binding$status=="PASS",hash(metadata_path)==binding$metadata_sha256,
 hash(file.path(root,"results/tables/bulk/frozen_multicellular_signature.tsv"))==binding$signature_sha256,
 hash(file.path(root,"scripts/export_retained_cell_counts.R"))==binding$exporter_sha256)
for(f in binding$files) {
 stopifnot(!grepl("[/\\\\]",f$path),file.info(file.path(input,f$path))$size==f$bytes,
  hash(file.path(input,f$path))==f$sha256)
}
md<-read_tsv(metadata_path)
stopifnot(nrow(md)==14846,!anyDuplicated(md$cell_id),all(md$passes_qc),all(md$scDblFinder.class=="singlet"),
 length(unique(md$sample_id))==4,length(unique(md$patient_id))==2,length(unique(md$cell_type))==9)
sig<-read_tsv(file.path(root,"results/tables/bulk/frozen_multicellular_signature.tsv"))
requested<-unique(c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5",
 "CKM","MYH1","MYH2","MYH7","TNNT3","TTN","H19","ADIPOQ","TENT5A","ESR1","METTL3","MSTN","SMAD3",sig$gene))
panels<-list(bulk_candidates=c("ANGPT1","LEP","ADIPOQ"),
 original_FAP_canonical_markers=c("PDGFRA","COL1A1","COL1A2","COL3A1","DCN","LUM","PI16","DPP4"),
 nominated_stromal_and_contamination_genes=c("CFD","APOD","POSTN","IGKC","ACTA1","MYL2","DES"))
biology<-unique(unlist(panels));union_targets<-unique(c(requested,biology))
features<-read_tsv(file.path(input,"cell_assay_features.tsv.gz"))
availability<-read_tsv(file.path(input,"cell_target_availability.tsv"))
assay_info<-read_tsv(file.path(input,"cell_assay_manifest.tsv"))
sparse<-read_tsv(file.path(input,"cell_target_counts.tsv.gz"))
stopifnot(length(requested)==158,length(biology)==18,length(union_targets)==165,
 !anyDuplicated(features[,c("assay","gene")]),!anyDuplicated(availability[,c("assay","gene")]),
 setequal(features$assay,c("RNA","SoupX")),setequal(assay_info$assay,c("RNA","SoupX")),nrow(assay_info)==2,
 nrow(availability)==length(union_targets)*2,setequal(availability$gene,union_targets),
 !anyDuplicated(sparse[,c("assay","cell_id","gene")]),all(sparse$assay%in%c("RNA","SoupX")),
 all(sparse$cell_id%in%md$cell_id),all(sparse$gene%in%union_targets),
 all(is.finite(sparse$UMI)),all(sparse$UMI>0),all(sparse$UMI==round(sparse$UMI)))
for(a in c("RNA","SoupX")) {
 f<-features$gene[features$assay==a];av<-availability[availability$assay==a,]
 stopifnot(nrow(av)==length(union_targets),setequal(av$gene,union_targets),
  identical(av$available,av$gene%in%f),identical(av$original_expression_requested,av$gene%in%requested),
  identical(av$closed_biology_target,av$gene%in%biology),
  length(f)==assay_info$n_genes[assay_info$assay==a],assay_info$n_cells[assay_info$assay==a]==nrow(md),
  all(sparse$gene[sparse$assay==a]%in%f),
  nrow(sparse[sparse$assay==a,])==assay_info$n_export_nonzero_records[assay_info$assay==a])
}
# All supplied genes are available in both actual assays. Unavailable features
# would be omitted from the intersection, never fabricated as zero counts.
original_availability<-do.call(rbind,lapply(c("RNA","SoupX"),function(a)data.frame(assay=a,gene=requested,
 available=requested%in%features$gene[features$assay==a])))
target<-Reduce(intersect,list(requested,features$gene[features$assay=="RNA"],features$gene[features$assay=="SoupX"]))
stopifnot(length(target)==158,all(biology%in%features$gene[features$assay=="RNA"]),all(biology%in%features$gene[features$assay=="SoupX"]))
dir.create(file.path(run,"analysis"),recursive=TRUE,showWarnings=FALSE)
run<-normalizePath(run,mustWork=TRUE);out<-file.path(run,"analysis")
dir.create(file.path(run,"results/tables/annotation"),recursive=TRUE)
stopifnot(file.copy(metadata_path,file.path(run,"results/tables/annotation/annotated_cell_metadata.tsv.gz")))
write_tsv<-function(x,n) {
 p<-file.path(out,n);con<-if(grepl("\\.gz$",p))gzfile(p,"wt") else file(p,"wt")
 on.exit(close(con));write.table(x,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
write_tsv(original_availability,"signature_requested_gene_availability.tsv")
write_tsv(do.call(rbind,lapply(names(panels),function(n)data.frame(panel=n,gene=panels[[n]]))),"signature_biology_target_panel_manifest.tsv")
libraries<-sort(unique(md$sample_id));lineages<-sort(unique(md$cell_type))
audit<-list();unexpected_sources<-list();assay_manifest<-list();bio_summaries<-list();bio_sources<-list()
unexpected<-intersect(c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5","ESR1"),target)
for(a in c("RNA","SoupX")) {
 message("Recompute full retained cell audit: ",a)
 z<-sparse[sparse$assay==a,]
 available<-union_targets[union_targets%in%features$gene[features$assay==a]]
 m<-sparseMatrix(i=match(z$gene,available),j=match(z$cell_id,md$cell_id),x=z$UMI,
  dims=c(length(available),nrow(md)),dimnames=list(available,md$cell_id))
 totals<-md[[paste0("nCount_",a)]]
 stopifnot(all(totals>0),all(totals==round(totals)),all(Matrix::colSums(m)<=totals))
 selected<-as.matrix(m[target,,drop=FALSE]);cp10k<-sweep(selected,2,pmax(totals,1),"/")*1e4;ln<-log1p(cp10k)
 assay_manifest[[a]]<-data.frame(assay=a,stage=if(a=="RNA")"before_ambient_correction" else "after_SoupX",
  n_genes=sum(features$assay==a),n_cells=nrow(md),n_target_genes=length(target),
  transformation="UMI_counts;per_cell_counts_divided_by_all_gene_UMIs_times_10000;natural_log1p",
  source_object="data/retained/cell_target_counts.tsv.gz")
 for(lib in libraries)for(ct in lineages) {
  ix<-which(md$sample_id==lib&md$cell_type==ct)
  if(!length(ix)) {
   audit[[paste(a,lib,ct)]]<-data.frame(assay=a,sample_id=lib,patient_id=sub("__.*$","",lib),side=sub("^.*__","",lib),
    cell_type=ct,gene=target,n_cells=0,n_detected=0,detection_fraction=NA,mean_UMI_per_cell=NA,mean_UMI_detected_cells=NA,
    mean_cp10k=NA,mean_log1p_cp10k=NA,median_log1p_cp10k=NA,p90_log1p_cp10k=NA,max_log1p_cp10k=NA,
    inferential_status="descriptive_library_lineage_distribution_no_cell_independence")
  } else {
   d<-selected[,ix,drop=FALSE];l<-ln[,ix,drop=FALSE];det<-rowSums(d>0)
   audit[[paste(a,lib,ct)]]<-data.frame(assay=a,sample_id=lib,patient_id=md$patient_id[ix[1]],side=md$side[ix[1]],
    cell_type=ct,gene=target,n_cells=length(ix),n_detected=det,detection_fraction=det/length(ix),mean_UMI_per_cell=rowMeans(d),
    mean_UMI_detected_cells=ifelse(det>0,rowSums(d)/det,NA),mean_cp10k=rowMeans(cp10k[,ix,drop=FALSE]),
    mean_log1p_cp10k=rowMeans(l),median_log1p_cp10k=apply(l,1,median),p90_log1p_cp10k=apply(l,1,quantile,.9),
    max_log1p_cp10k=apply(l,1,max),inferential_status="descriptive_library_lineage_distribution_no_cell_independence")
  }
 }
 unexpected_sources[[a]]<-do.call(rbind,lapply(unexpected,function(g)data.frame(assay=a,cell_id=md$cell_id,sample_id=md$sample_id,
  patient_id=md$patient_id,side=md$side,cell_type=md$cell_type,gene=g,UMI=selected[g,],log1p_cp10k=ln[g,])))
 b<-as.matrix(m[biology,,drop=FALSE]);bn<-sweep(b,2,pmax(totals,1),"/")*1e4;bl<-log1p(bn)
 bio_sources[[a]]<-do.call(rbind,lapply(biology,function(g)data.frame(assay=a,cell_id=md$cell_id,sample_id=md$sample_id,
  patient_id=md$patient_id,side=md$side,cell_type=md$cell_type,gene=g,UMI=b[g,],all_gene_total_UMI=totals,cp10k=bn[g,],log1p_cp10k=bl[g,])))
 for(lib in libraries)for(ct in lineages)for(g in biology) {
  ix<-which(md$sample_id==lib&md$cell_type==ct);v<-b[g,ix];n<-length(ix);total<-sum(v);s<-sort(v,decreasing=TRUE)
  bio_summaries[[paste(a,lib,ct,g)]]<-data.frame(assay=a,sample_id=lib,patient_id=sub("__.*$","",lib),side=sub("^.*__","",lib),cell_type=ct,gene=g,
   n_cells=n,feature_available=TRUE,n_detected=sum(v>0),total_UMI=total,detection_fraction=if(n)mean(v>0) else NA,
   mean_UMI_per_cell=if(n)mean(v) else NA,mean_cp10k=if(n)mean(bn[g,ix]) else NA,mean_log1p_cp10k=if(n)mean(bl[g,ix]) else NA,
   median_log1p_cp10k=if(n)median(bl[g,ix]) else NA,fraction_UMI_in_top1_cell=if(total>0)s[1]/total else NA,
   fraction_UMI_in_top5_cells=if(total>0)sum(head(s,5))/total else NA,
   fraction_UMI_in_top1percent_cells=if(total>0)sum(head(s,max(1,ceiling(n*.01))))/total else NA)
 }
}
audit<-do.call(rbind,audit);expression<-do.call(rbind,unexpected_sources);bio<-do.call(rbind,bio_summaries);bs<-do.call(rbind,bio_sources)
write_tsv(audit,"signature_expression_by_library_lineage_before_after.tsv")
write_tsv(do.call(rbind,assay_manifest),"signature_expression_assay_manifest.tsv")
write_tsv(expression,"signature_unexpected_gene_cell_source.tsv.gz")
key<-c("sample_id","patient_id","side","cell_type","gene","n_cells")
before<-audit[audit$assay=="RNA",];after<-audit[audit$assay=="SoupX",]
comparison<-merge(before[,c(key,"detection_fraction","mean_UMI_per_cell","mean_log1p_cp10k")],
 after[,c(key,"detection_fraction","mean_UMI_per_cell","mean_log1p_cp10k")],by=key,suffixes=c(".before",".after"))
comparison$detection_fraction_change<-comparison$detection_fraction.after-comparison$detection_fraction.before
comparison$mean_UMI_change<-comparison$mean_UMI_per_cell.after-comparison$mean_UMI_per_cell.before
write_tsv(comparison,"signature_expression_correction_changes.tsv")
write_tsv(merge(sig,audit[audit$assay=="SoupX",],by=c("gene","cell_type"),all.x=TRUE),"signature_frozen_gene_expression_all_four_libraries.tsv")
fields<-c("cell_id","nCount_RNA","nFeature_RNA","nCount_SoupX","nFeature_SoupX","subsets_Mito_percent","subsets_Ribo_percent",
 "scDblFinder.score","scDblFinder.class","passes_qc_metrics","passes_qc","seurat_clusters","annotation_confidence")
positive<-expression[expression$UMI>0,];ix<-match(positive$cell_id,md$cell_id)
trace<-cbind(positive,md[ix,setdiff(fields,"cell_id"),drop=FALSE]);write_tsv(trace,"signature_positive_cell_QC_trace.tsv.gz")
keys<-c("assay","sample_id","patient_id","side","cell_type","gene");groups<-unique(expression[,keys])
qc_summaries<-lapply(seq_len(nrow(groups)),function(i) {
 g<-groups[i,];keep<-rep(TRUE,nrow(expression));for(k in keys)keep<-keep&expression[[k]]==g[[k]]
 d<-expression[keep,];p<-d[d$UMI>0,];info<-md[match(p$cell_id,md$cell_id),];total<-sum(d$UMI)
 cbind(g,data.frame(n_cells=nrow(d),n_detected=nrow(p),total_gene_UMI=total,largest_single_cell_UMI=if(nrow(p))max(p$UMI) else 0,
  largest_cell_fraction_of_gene_UMIs=if(total>0)max(p$UMI)/total else NA_real_,
  positive_cells_median_nCount_RNA=if(nrow(p))median(info$nCount_RNA) else NA_real_,
  positive_cells_max_mito_percent=if(nrow(p))max(info$subsets_Mito_percent) else NA_real_,
  positive_cells_max_original_doublet_score=if(nrow(p))max(info$scDblFinder.score) else NA_real_,
  all_positive_cells_original_QC_singlets=if(nrow(p))all(info$passes_qc&info$scDblFinder.class=="singlet") else NA,
  interpretation="descriptive_original_QC_status_does_not_authenticate_lineage_transcription"))
})
write_tsv(do.call(rbind,qc_summaries),"signature_positive_cell_QC_concentration.tsv")
bio$scientific_scope<-"descriptive_two_donor_closed_target_audit_no_cell_level_independence_no_signature_training"
write_tsv(bio,"signature_biology_target_gene_expression_all_libraries_lineages.tsv")
write_tsv(bs,"signature_biology_all_target_cell_source.tsv.gz")
write_tsv(bs[bs$gene%in%panels$bulk_candidates,],"signature_biology_bulk_candidate_cell_source.tsv.gz")
shares<-do.call(rbind,lapply(split(bio,interaction(bio$assay,bio$sample_id,bio$gene,drop=TRUE)),function(z) {
 z$fraction_gene_UMI_in_lineage<-if(sum(z$total_UMI,na.rm=TRUE)>0)z$total_UMI/sum(z$total_UMI,na.rm=TRUE) else NA_real_
 z$fraction_detected_cells_in_lineage<-if(sum(z$n_detected,na.rm=TRUE)>0)z$n_detected/sum(z$n_detected,na.rm=TRUE) else NA_real_;z
}))
write_tsv(shares[,c("assay","sample_id","cell_type","gene","n_cells","n_detected","total_UMI","fraction_gene_UMI_in_lineage","fraction_detected_cells_in_lineage")],"signature_biology_retained_cell_expression_lineage_provenance.tsv")
pairs<-do.call(rbind,lapply(split(bio,interaction(bio$assay,bio$cell_type,bio$gene,bio$patient_id,drop=TRUE)),function(z) {
 cc<-z[z$side=="concave",];cv<-z[z$side=="convex",];stopifnot(nrow(cc)==1,nrow(cv)==1)
 data.frame(assay=cc$assay,cell_type=cc$cell_type,gene=cc$gene,patient_id=cc$patient_id,n_concave=cc$n_cells,n_convex=cv$n_cells,
  concave_mean_UMI=cc$mean_UMI_per_cell,convex_mean_UMI=cv$mean_UMI_per_cell,mean_UMI_difference=cc$mean_UMI_per_cell-cv$mean_UMI_per_cell,
  mean_log1p_cp10k_difference=cc$mean_log1p_cp10k-cv$mean_log1p_cp10k,detection_fraction_difference=cc$detection_fraction-cv$detection_fraction)
}))
write_tsv(pairs,"signature_biology_target_gene_descriptive_within_donor_differences.tsv")
manifest<-data.frame(field=c("input_type","input_relative_path","input_sha256","metadata_sha256","rebuild_script_sha256",
 "R_version","n_retained_singlets","n_libraries","n_lineages","n_requested_genes","n_available_genes_in_both_assays","n_summary_rows","scientific_scope"),
 value=c("genuine_retained_sparse_target_UMIs_with_full_assay_totals","data/retained/cell_target_counts.tsv.gz",hash(file.path(input,"cell_target_counts.tsv.gz")),
 hash(metadata_path),hash(script),as.character(getRversion()),nrow(md),length(libraries),length(lineages),length(requested),length(target),nrow(audit),
 "fresh_sparse_count_reconstruction;no_new_Seurat_or_RDS_reconstruction;no_cell_level_hypothesis_test"))
write_tsv(manifest,"signature_cell_expression_run_manifest.tsv")
write_tsv(data.frame(field=c("input_type","input_sha256","n_retained_cells","n_gene_targets","n_summary_rows","n_candidate_cell_source_rows","n_all_target_cell_source_rows","R_version","scope"),
 value=c("genuine_retained_sparse_target_UMIs",hash(file.path(input,"cell_target_counts.tsv.gz")),nrow(md),length(biology),nrow(bio),sum(bs$gene%in%panels$bulk_candidates),nrow(bs),as.character(getRversion()),
 "fresh_sparse_count_reconstruction;two_donors_descriptive_only;historical_RDS_receipt_not_reused")),"signature_biology_execution_manifest.tsv")
report<-list(status="COMPLETED_PENDING_INDEPENDENT_COMPARISON",scope="All158 requested-gene summaries plus fixed18 biology counts reconstructed from genuine sparse cell inputs; no old summary is read",
 input_files=binding$files,metadata_sha256=hash(metadata_path),rebuild_script_sha256=hash(script),
 n_cells=nrow(md),n_union_targets=length(union_targets),n_nonzero_count_records=nrow(sparse),
 n_expression_summaries=nrow(audit),n_unexpected_cell_records=nrow(expression),n_positive_QC_records=nrow(trace),
 n_QC_concentration_groups=length(qc_summaries),n_biology_summary_groups=nrow(bio),n_biology_cell_records=nrow(bs),
 n_biology_donor_contrasts=nrow(pairs),n_biology_share_groups=nrow(shares),
 historical_export_source_object_sha256=binding$input_object_sha256,
 historical_export_scope="Object source of supplied count export only; this run does not load or newly reconstruct that object")
jsonlite::write_json(report,file.path(run,"cell_audit_replay.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat(run,"\n",sep="")
