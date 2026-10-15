#!/usr/bin/env Rscript
# Available corrected-count fallback when external per-cell RDS files cannot read.
# These summaries do not estimate per-cell detection frequency.
out<-"analysis"
read_tsv<-function(p) read.delim(if(grepl("\\.gz$",p))gzfile(p) else p,check.names=FALSE)
write_tsv<-function(x,n) write.table(x,file.path(out,paste0("signature_",n,".tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
sig<-read_tsv("results/tables/bulk/frozen_multicellular_signature.tsv")
cells<-read_tsv("results/tables/pseudobulk/cell_counts_per_pseudobulk.tsv")
de<-read_tsv("results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz")
targets<-unique(c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5",
 "H19","ADIPOQ","TENT5A","ESR1","METTL3","MSTN","SMAD3",sig$gene))
rows<-list()
for(ct in unique(de$cell_type)) {
 file<-file.path("results/tables/pseudobulk/counts",paste0(gsub("[^A-Za-z0-9]+","_",ct),".tsv.gz"))
 tab<-read_tsv(file);m<-as.matrix(tab[,-1,drop=FALSE]);rownames(m)<-tab$gene
 totals<-colSums(m)
 for(lib in colnames(m))for(g in targets) {
  n<-cells$n_cells[cells$cell_type==ct&cells$pb_id==lib]
  exists<-g%in%rownames(m)
  v<-if(exists)m[g,lib] else NA_real_
  rows[[paste(ct,lib,g)]]<-data.frame(gene=g,cell_type=ct,sample_id=lib,patient_id=sub("__.*$","",lib),side=sub("^.*__","",lib),
   n_singlets=n,feature_present=exists,total_corrected_UMI=v,mean_corrected_UMI_per_cell=v/n,
   pseudobulk_total_all_gene_UMI=totals[lib],corrected_pseudobulk_CPM=v/totals[lib]*1e6,
   per_cell_detection_fraction=NA_real_,in_lineage_frozen_signature=any(sig$gene==g&sig$cell_type==ct),
   expression_source="retained_SoupX_corrected_pseudobulk_count_sum",
   detection_status="cannot_recover_detection_frequency_from_pseudobulk_sums")
 }
}
rows<-do.call(rbind,rows)
write_tsv(rows,"corrected_pseudobulk_gene_expression_four_libraries")
# Pooled corrected cell-level summaries exist locally for the original marker panel.
# Restrict to requested DES; pooling across libraries is labelled and cannot show
# the within-library before/after distributions the reviewers requested.
marker<-read_tsv("results/tables/annotation/marker_specificity_by_cell_type.tsv")
marker<-marker[marker$gene%in%c("IGKC","ACTA1","MYL2","DES","CFD","APOD"),]
marker$sample_scope<-"all_four_discovery_libraries_pooled"
marker$assay<-"SoupX_corrected"
marker$interpretation<-"descriptive_pooled_marker_summary_not_library_specific_no_pre_post_comparison"
write_tsv(marker,"unexpected_genes_retained_pooled_cell_summary")
cell_status<-"not_completed_in_this_output_directory"
cell_details<-"Run signature_cell_expression_audit.R on an annotated object; no detection values inferred from count sums"
cell_files<-file.path(out,c("signature_cell_expression_run_manifest.tsv",
 "signature_expression_by_library_lineage_before_after.tsv"))
if(all(file.exists(cell_files))) {
 manifest<-read_tsv(cell_files[1]);values<-setNames(manifest$value,manifest$field)
 actual<-read_tsv(cell_files[2])
 stopifnot(values[["n_retained_singlets"]]=="14846",values[["n_libraries"]]=="4",
           values[["n_lineages"]]=="9",nrow(actual)==as.integer(values[["n_summary_rows"]]),
           setequal(actual$assay,c("RNA","SoupX")),length(unique(actual$sample_id))==4,
           length(unique(actual$cell_type))==9)
 cell_status<-"completed_real_annotated_singlets"
 cell_details<-paste0(nrow(actual),"_real_gene_lineage_library_assay_summaries; ",
                    "source_object_MD5=",values[["input_object_md5"]],
                    "; descriptive_only; does_not_authenticate_lineage_biology")
}
write_tsv(data.frame(request=c("per_cell_before_after_all_four_libraries_all_nine_lineages","corrected_pseudobulk_mean_expression_all_four_libraries_five_eligible_lineages","pooled_corrected_DES_detection_all_nine_lineages"),
 status=c(cell_status,"completed","completed"),
 details=c(cell_details,
 paste0(length(targets),"_gene_targets_with_feature_availability_and_no_detection_inference"),"only_DES_is_present_in_retained_canonical_marker_panel")),"cell_expression_availability")
message("Available corrected-expression summaries written without fabricated per-cell detection values.")
