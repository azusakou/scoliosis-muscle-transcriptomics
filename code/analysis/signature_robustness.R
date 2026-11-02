#!/usr/bin/env Rscript
source("scripts/analysis_paths.R")
suppressPackageStartupMessages(library(stats))
out <- "analysis"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_tsv <- function(p) {p <- analysis_input(p); read.delim(if(grepl("\\.gz$",p)) gzfile(p) else p, check.names=FALSE)}
write_tsv <- function(x,n) write.table(x,file.path(out,paste0("signature_",n,".tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
read_matrix <- function(p) { x<-read_tsv(p); m<-as.matrix(x[,-1,drop=FALSE]); rownames(m)<-x[[1]]; storage.mode(m)<-"double"; m }
sig <- read_tsv("results/tables/bulk/frozen_multicellular_signature.tsv")
conflicts <- unique(read_tsv("results/tables/bulk/frozen_multicellular_signature_direction_conflicts.tsv")$gene)
source_sig <- sig
sig <- sig[!sig$gene %in% conflicts,,drop=FALSE]
cohorts <- list(
 GSE254300_AIS=list(x=read_matrix("results/tables/bulk/GSE254300_AIS_label_blind_projection_logCPM.tsv.gz"),
                    info=read_tsv("results/tables/bulk/multicellular_axis_scores.tsv")),
 GSE311180_AIS=list(x=read_matrix("results/tables/bulk_gse311180_requant/GSE311180_label_blind_projection_logCPM.tsv.gz"),
                    info=read_tsv("results/tables/bulk_gse311180_requant/frozen_axis_scores.tsv")))
for(k in names(cohorts)) {
 x<-cohorts[[k]]$x; z<-t(scale(t(x))); z[!is.finite(z)]<-0
 cohorts[[k]]$z<-z
}
exact_test <- function(d) {
 obs<-mean(d); signs<-as.matrix(expand.grid(rep(list(c(-1,1)),length(d))))
 mean(abs(as.numeric(signs%*%d/length(d)))>=abs(obs)-1e-12)
}
paired <- function(score,info) {
 ids<-sort(unique(info$patient_id)); vapply(ids,function(id) {
  cc<-info$sample_id[info$patient_id==id & info$side=="concave"]
  cv<-info$sample_id[info$patient_id==id & info$side=="convex"]
  stopifnot(length(cc)==1,length(cv)==1); as.numeric(score[cc]-score[cv])
 },numeric(1))
}
score <- function(s,z) {
 w<-aggregate(weight~gene,s,mean); w<-w[is.finite(w$weight)&w$weight!=0&w$gene%in%rownames(z),,drop=FALSE]
 if(!nrow(w)) return(list(score=setNames(rep(NA_real_,ncol(z)),colnames(z)),n_genes=0,weights=w))
 list(score=colSums(z[w$gene,,drop=FALSE]*w$weight)/sum(abs(w$weight)),n_genes=nrow(w),weights=w)
}
summarise <- function(d,label,cohort,n_genes) {
 data.frame(analysis=label,cohort=cohort,n_pairs=length(d),n_projected_unique_genes=n_genes,
  mean_paired_difference=mean(d),sd_paired_difference=sd(d),positive_pairs=sum(d>0),
  exact_two_sided_sign_flip_p=exact_test(d))
}

# Broad, explicit gene exclusions. No exclusion list is fitted to bulk side effects.
deconv <- read_tsv("results/tables/bulk/deconvolution_marker_reference.tsv")$gene
canon_lines <- readLines("config/canonical_markers.yaml")
canon <- unique(unlist(strsplit(gsub(".*\\[|\\].*","",canon_lines),",[ ]*")))
top_markers <- read_tsv("results/tables/annotation/top15_cluster_markers.tsv")$gene
all_markers <- read_tsv("results/tables/annotation/cluster_differential_markers.tsv.gz")
# Annotation tests only define an exclusion list; their cell-level P values are not
# evidence for side-associated biology. This sensitivity deliberately removes all
# positively enriched annotation markers, regardless of their nominal P values.
all_markers <- unique(all_markers$gene[all_markers$avg_log2FC>0])
suspect <- unique(c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5",
 "CKM","MYH1","MYH2","MYH7","TNNT1","TNNT3","TPM2","TTN","ACTC1",
 "IGLC1","IGLC2","IGLC3","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2","IGHM","JCHAIN"))
sets <- list(contamination_prone=suspect,deconvolution_markers=unique(deconv),
 canonical_and_deconvolution_markers=unique(c(canon,deconv)),
 all_annotation_and_deconvolution_markers=unique(c(all_markers,canon,deconv)),
 contamination_and_all_markers=unique(c(suspect,all_markers,canon,deconv)))
exclude_manifest <- do.call(rbind,lapply(names(sets),function(n) data.frame(exclusion_set=n,gene=sort(sets[[n]]),in_original_signature=sort(sets[[n]])%in%sig$gene)))
write_tsv(exclude_manifest,"gene_exclusion_manifest")
variants<-list(original=sig)
for(ct in sort(unique(sig$cell_type))) variants[[paste0("leave_out__",ct)]]<-sig[sig$cell_type!=ct,,drop=FALSE]
for(ct in sort(unique(sig$cell_type))) variants[[paste0("lineage_only__",ct)]]<-sig[sig$cell_type==ct,,drop=FALSE]
for(n in names(sets)) variants[[paste0("remove__",n)]]<-sig[!sig$gene%in%sets[[n]],,drop=FALSE]
cell_counts<-read_tsv("results/tables/pseudobulk/cell_counts_per_pseudobulk.tsv")
eligibility<-do.call(rbind,lapply(split(cell_counts,cell_counts$cell_type),function(x) data.frame(cell_type=x$cell_type[1],
 min_cells_per_library=min(x$n_cells),max_cells_per_library=max(x$n_cells),eligible_at_20=all(x$n_cells>=20),eligible_at_50=all(x$n_cells>=50),eligible_at_100=all(x$n_cells>=100))))
write_tsv(eligibility,"pseudobulk_cell_threshold_eligibility")
variants$minimum_50_cells_all_libraries<-sig[sig$cell_type%in%eligibility$cell_type[eligibility$eligible_at_50],,drop=FALSE]
variants$minimum_100_cells_all_libraries<-sig[sig$cell_type%in%eligibility$cell_type[eligibility$eligible_at_100],,drop=FALSE]
summary_all<-list(); effects_all<-list(); score_all<-list(); usage_all<-list()
for(v in names(variants)) {
 d_pool<-numeric(); n_pool<-character()
 for(k in names(cohorts)) {
  pr<-score(variants[[v]],cohorts[[k]]$z)
  if(pr$n_genes==0) next
  d<-paired(pr$score,cohorts[[k]]$info)
  summary_all[[paste(v,k)]]<-summarise(d,v,k,as.character(pr$n_genes))
  effects_all[[paste(v,k)]]<-data.frame(analysis=v,cohort=k,patient_id=names(d),paired_difference=as.numeric(d))
  info<-cohorts[[k]]$info; info$score<-pr$score[info$sample_id]; info$analysis<-v; info$cohort<-k
  score_all[[paste(v,k)]]<-info[,c("analysis","cohort","sample_id","patient_id","side","score")]
  usage_all[[paste(v,k)]]<-data.frame(analysis=v,cohort=k,pr$weights)
  d_pool<-c(d_pool,d);n_pool<-c(n_pool,paste0(k,":",pr$n_genes))
 }
 if(length(d_pool)==12) summary_all[[paste(v,"pooled")]]<-summarise(d_pool,v,"pooled_12_AIS",paste(n_pool,collapse=";"))
}
summary_all<-do.call(rbind,summary_all)
summary_all$inferential_status<-"post_hoc_sensitivity_no_primary_reselection"
summary_all$original_20_gene_projection_gate_met<-vapply(strsplit(summary_all$n_projected_unique_genes,";",fixed=TRUE),function(x) all(as.numeric(sub(".*:","",x))>=20),logical(1))
summary_all$FDR_five_lineage_subscores<-NA_real_
for(k in unique(summary_all$cohort)) {
 ix<-grepl("^lineage_only__",summary_all$analysis)&summary_all$cohort==k
 summary_all$FDR_five_lineage_subscores[ix]<-p.adjust(summary_all$exact_two_sided_sign_flip_p[ix],"BH")
}
write_tsv(summary_all,"bulk_sensitivity_tests")
effects_all<-do.call(rbind,effects_all);write_tsv(effects_all,"bulk_sensitivity_patient_effects")
write_tsv(do.call(rbind,score_all),"bulk_sensitivity_sample_scores")
write_tsv(do.call(rbind,usage_all),"bulk_sensitivity_projected_weights")
lineage_cor<-list()
for(k in names(cohorts)) {
 e<-effects_all[effects_all$cohort==k&grepl("^lineage_only__",effects_all$analysis),]
 wide<-reshape(e[,c("patient_id","analysis","paired_difference")],idvar="patient_id",timevar="analysis",direction="wide")
 m<-as.matrix(wide[,-1,drop=FALSE]);colnames(m)<-sub("paired_difference.lineage_only__","",colnames(m),fixed=TRUE)
 pairs<-combn(colnames(m),2,simplify=FALSE)
 lineage_cor[[k]]<-do.call(rbind,lapply(pairs,function(p) data.frame(cohort=k,lineage_1=p[1],lineage_2=p[2],n_pairs=nrow(m),spearman_r=cor(m[,p[1]],m[,p[2]],method="spearman"),interpretation="descriptive_correlated_gene_subscores_not_coordination_test")))
}
write_tsv(do.call(rbind,lineage_cor),"lineage_subscore_correlations")

# Discovery null: complete paired side flips, shared across all genes/lineages.
de<-read_tsv("results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz")
effect_cols<-grep("^paired_logCPM_difference",names(de),value=TRUE)
stopifnot(length(effect_cols)==2)
flip_grid<-as.matrix(expand.grid(H01045=c(1,-1),H66206=c(1,-1)))
concordance<-do.call(rbind,lapply(split(de,de$cell_type),function(x) {
 a<-as.matrix(x[,effect_cols]);nz<-rowSums(abs(a)>1e-8)==2;n<-sum(nz);cct<-sum(a[nz,1]*a[nz,2]>0)
 data.frame(cell_type=x$cell_type[1],n_tested_genes=nrow(x),n_nonzero_two_donor_genes=n,
 observed_concordant=cct,observed_fraction=cct/n,independent_sign_null_expected_count=n/2,
 numerical_zero_tolerance=1e-8,original_exact_sign_concordant=sum(x$donor_sign_concordant),
 paired_flip_null_min_count=min(cct,n-cct),paired_flip_null_max_count=max(cct,n-cct),
 interpretation="50_percent_descriptive_expectation_gene_dependence_precludes_binomial_inference")
}))
write_tsv(concordance,"discovery_concordance_counts")
pipeline_signatures<-list();pipeline_tests<-list();pipeline_effects<-list();pipeline_counts<-list()
for(i in seq_len(nrow(flip_grid))) {
 d<-de;a<-sweep(as.matrix(d[,effect_cols]),2,flip_grid[i,],"*")
 d$mean_donor_effect<-rowMeans(a);d$donor_sign_concordant<-a[,1]*a[,2]>0
 d<-d[d$donor_sign_concordant&!grepl("^(MT-|RPS|RPL|HBA|HBB)",d$gene),,drop=FALSE]
 s<-do.call(rbind,lapply(split(d,d$cell_type),function(z) {z<-z[order(abs(z$mean_donor_effect),decreasing=TRUE),];z<-head(z,30);z$rank_within_lineage<-seq_len(nrow(z));z}))
 s$weight<-sign(s$mean_donor_effect);s$configuration<-paste(flip_grid[i,],collapse=",")
 pipeline_signatures[[i]]<-s[,c("configuration","gene","cell_type","rank_within_lineage","mean_donor_effect","weight")]
 w<-aggregate(weight~gene,s,mean)
 pipeline_counts[[i]]<-data.frame(configuration=s$configuration[1],n_lineage_entries=nrow(s),n_unique_genes=nrow(w),n_zero_weight_conflicts=sum(w$weight==0),n_usable_genes=sum(w$weight!=0))
 dp<-numeric()
 for(k in names(cohorts)) {
  pr<-score(s,cohorts[[k]]$z);d<-paired(pr$score,cohorts[[k]]$info)
  pipeline_tests[[paste(i,k)]]<-summarise(d,s$configuration[1],k,as.character(pr$n_genes))
  pipeline_effects[[paste(i,k)]]<-data.frame(configuration=s$configuration[1],cohort=k,patient_id=names(d),paired_difference=d)
  dp<-c(dp,d)
 }
 pipeline_tests[[paste(i,"pool")]]<-summarise(dp,s$configuration[1],"pooled_12_AIS","cohort_specific")
}
pipeline_tests<-do.call(rbind,pipeline_tests)
null_summary<-do.call(rbind,lapply(split(pipeline_tests,pipeline_tests$cohort),function(z) {
 obs<-z$mean_paired_difference[z$analysis=="1,1"]
 data.frame(cohort=z$cohort[1],n_paired_discovery_configurations=nrow(z),observed_mean=obs,
 exact_pipeline_two_sided_tail=mean(abs(z$mean_paired_difference)>=abs(obs)-1e-12),
 minimum_possible_two_sided_tail=.5,
 interpretation="four_global_pair_flips_two_unsigned_selections_no_cell_or_gene_independence")
}))
write_tsv(do.call(rbind,pipeline_signatures),"pipeline_flip_reselected_signatures")
write_tsv(do.call(rbind,pipeline_counts),"pipeline_flip_signature_counts")
write_tsv(pipeline_tests,"pipeline_flip_projection_tests")
write_tsv(do.call(rbind,pipeline_effects),"pipeline_flip_patient_effects")
write_tsv(null_summary,"pipeline_flip_null_summary")
selected_original<-pipeline_signatures[[1]]
selection_check<-merge(source_sig[,c("gene","cell_type","rank_within_lineage","weight","mean_donor_effect")],
 selected_original[,c("gene","cell_type","rank_within_lineage","weight","mean_donor_effect")],
 by=c("gene","cell_type"),suffixes=c(".retained",".reselected"),all=TRUE)
stopifnot(nrow(selection_check)==150,all(complete.cases(selection_check)),
 all(selection_check$rank_within_lineage.retained==selection_check$rank_within_lineage.reselected),
 all(selection_check$weight.retained==selection_check$weight.reselected),
 max(abs(selection_check$mean_donor_effect.retained-selection_check$mean_donor_effect.reselected))<1e-10)
write_tsv(data.frame(n_lineage_gene_entries=nrow(selection_check),identical_gene_lineage_membership=TRUE,
 identical_within_lineage_ranks=TRUE,identical_sign_weights=TRUE,
 max_abs_mean_effect_difference=max(abs(selection_check$mean_donor_effect.retained-selection_check$mean_donor_effect.reselected))),
 "pipeline_original_selection_reconstruction_check")

# Literature-nominated genes: report actual availability and both donor effects.
genes<-c("H19","ADIPOQ","TENT5A","ESR1","METTL3","MSTN","SMAD3")
literature_sc<-do.call(rbind,lapply(sort(unique(de$cell_type)),function(ct) do.call(rbind,lapply(genes,function(g) {
 z<-de[de$cell_type==ct&de$gene==g,,drop=FALSE]
 if(!nrow(z)) return(data.frame(gene=g,cell_type=ct,tested=FALSE,edgeR_logFC=NA,logCPM=NA,donor_H01045=NA,donor_H66206=NA,concordant_original_exact_sign=NA,both_donor_effects_above_1e_minus8=NA,concordant_tolerance_1e_minus8=NA,in_frozen_signature=any(source_sig$gene==g&source_sig$cell_type==ct)))
 a<-as.numeric(z[,effect_cols]);larger<-all(abs(a)>1e-8)
 data.frame(gene=g,cell_type=ct,tested=TRUE,edgeR_logFC=z$logFC,logCPM=z$logCPM,donor_H01045=z[[effect_cols[1]]],donor_H66206=z[[effect_cols[2]]],concordant_original_exact_sign=z$donor_sign_concordant,both_donor_effects_above_1e_minus8=larger,concordant_tolerance_1e_minus8=larger&&(all(a>0)||all(a<0)),in_frozen_signature=any(source_sig$gene==g&source_sig$cell_type==ct))
}))))
write_tsv(literature_sc,"literature_nominated_discovery_effects")
bulk_de<-list(GSE254300_AIS="results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz",
 GSE254300_CS="results/tables/disease_control_gse254300_cs/GSE254300_CS_paired_edgeR_full.tsv.gz",
 GSE311180_AIS="results/tables/bulk_gse311180_requant/GSE311180_paired_edgeR_full.tsv.gz")
bulk_gene<-do.call(rbind,lapply(names(bulk_de),function(k) {
 tab<-read_tsv(bulk_de[[k]]);do.call(rbind,lapply(genes,function(g) {
  data_gene<-if(g=="TENT5A"&&!g%in%tab$gene&&"FAM46A"%in%tab$gene)"FAM46A" else g
  z<-tab[tab$gene==data_gene,,drop=FALSE]
  data.frame(cohort=k,gene=g,input_gene_symbol=data_gene,tested=nrow(z)>0,logFC=if(nrow(z))z$logFC else NA,
   logCPM=if(nrow(z))z$logCPM else NA,PValue=if(nrow(z))z$PValue else NA,FDR_original_gene_family=if(nrow(z))z$FDR else NA,
   interpretation="targeted_post_hoc_report_original_genomewide_FDR_no_new_claim_of_replication")
 }))
}))
write_tsv(bulk_gene,"literature_nominated_bulk_effects")

# Cross-check full reconstruction against retained original patient-score output.
checks<-do.call(rbind,lapply(names(cohorts),function(k) {
 pr<-score(sig,cohorts[[k]]$z);old<-cohorts[[k]]$info
 data.frame(cohort=k,n_projected_genes=pr$n_genes,max_abs_sample_score_difference=max(abs(pr$score[old$sample_id]-old$score)))
}))
stopifnot(all(checks$max_abs_sample_score_difference<1e-10))
write_tsv(checks,"original_score_reconstruction_check")
message("Signature robustness outputs written; original scores reproduced.")
