#!/usr/bin/env Rscript
# Refit only the original CS paired and AIS-by-side interaction gene models.
specificity_hash <- function(p)digest::digest(file=p,algo="sha256")
specificity_read <- function(p)read.delim(if(grepl("[.]gz$",p))gzfile(p) else p,check.names=FALSE,stringsAsFactors=FALSE)
specificity_matrix <- function(p) {
  x<-specificity_read(p)
  if(!nrow(x)||anyDuplicated(x[[1]])||anyDuplicated(names(x))||any(!nzchar(x[[1]])))stop("Invalid specificity count identifiers")
  m<-as.matrix(x[,-1,drop=FALSE]);rownames(m)<-x[[1]];storage.mode(m)<-"double"
  if(any(!is.finite(m))||any(m<0)||any(m!=round(m)))stop("Invalid specificity count values")
  m
}
specificity_close <- function(x,y,label,tolerance=1e-7) {
 if(length(x)!=length(y)||any(!is.finite(x))||any(!is.finite(y))||any(abs(x-y)>tolerance*pmax(1,abs(y))))stop("Specificity model mismatch: ",label)
}
specificity_reconstruct <- function(A,extra,CS,genes,info) {
  ais<-info$cohort_label=="AIS";cs<-info$cohort_label=="CS"
  if(anyDuplicated(genes)||length(genes)!=14967L||nrow(A)!=14296L||nrow(extra)!=671L||
     !identical(rownames(CS),genes)||!setequal(c(rownames(A),rownames(extra)),genes)||
     length(intersect(rownames(A),rownames(extra)))||
     !identical(colnames(A),info$processed_sample_id[ais])||!identical(colnames(extra),colnames(A))||
     !identical(colnames(CS),info$processed_sample_id[cs]))stop("Changed specificity gene/sample partition")
  m<-matrix(NA_real_,length(genes),nrow(info),dimnames=list(genes,info$processed_sample_id))
  m[rownames(A),ais]<-A;m[rownames(extra),ais]<-extra;m[,cs]<-CS
  if(any(!is.finite(m))||any(m<0)||any(m!=round(m)))stop("Incomplete specificity joint counts")
  m
}
specificity_compare <- function(actual,expected,label) {
  if(!identical(actual$gene,expected$gene)||!identical(names(actual),names(expected)))stop("Changed specificity model gene/schema: ",label)
  differences<-list()
  for(field in c("logFC","logCPM","F","PValue","FDR")) {
    specificity_close(actual[[field]],expected[[field]],paste(label,field))
    if(field%in%c("PValue","FDR") && (!all(actual[[field]]>=0 & actual[[field]]<=1)||
       !identical(actual[[field]]<.05,expected[[field]]<.05)))stop("Changed probability decision: ",label,"/",field)
    if(field=="logFC"&&!identical(sign(actual[[field]]),sign(expected[[field]])))stop("Changed gene effect sign: ",label)
    differences[[field]]<-list(max_abs=max(abs(actual[[field]]-expected[[field]])),
      max_scaled=max(abs(actual[[field]]-expected[[field]])/pmax(1,abs(expected[[field]]))))
  }
  for(field in setdiff(names(actual),c("logFC","logCPM","F","PValue","FDR")))
    if(!identical(actual[[field]],expected[[field]]))stop("Changed specificity identifier/count: ",field)
  differences
}
specificity_main <- function() {
 suppressPackageStartupMessages(library(edgeR))
 args<-commandArgs(trailingOnly=TRUE);allowed<-c("--root","--input-dir","--output-dir");flags<-sub("=.*$","",args)
 if(any(!flags%in%allowed)||any(!grepl("=.+$",args))||anyDuplicated(flags))stop("Unknown, duplicate or empty CLI flag")
 value<-function(flag,default){z<-args[startsWith(args,paste0(flag,"="))];if(length(z))sub(paste0("^",flag,"="),"",z)else default}
 script<-normalizePath(sub("^--file=","",commandArgs(FALSE)[grep("^--file=",commandArgs(FALSE))][1]),mustWork=TRUE)
 root<-normalizePath(value("--root",dirname(dirname(script))),mustWork=TRUE)
 input<-normalizePath(value("--input-dir",file.path(root,"data/retained/gse254300_specificity")),mustWork=TRUE)
 logs<-file.path(root,"logs/specificity_model_replay");dir.create(logs,recursive=TRUE,showWarnings=FALSE)
 run<-value("--output-dir",tempfile("run-",tmpdir=logs));if(file.exists(run))stop("Specificity output already exists; use a new directory")
 binding_path<-file.path(input,"specificity_counts_export.json");binding<-jsonlite::fromJSON(binding_path,simplifyVector=FALSE)
 md_path<-file.path(root,"metadata/audit/bulk_samples.tsv")
 A_path<-file.path(root,"analysis/bulk_gate_inputs/GSE254300_tested_counts.tsv.gz")
 producer_path<-file.path(root,"scripts/61_gse254300_cs_specificity.R")
 cs_path<-file.path(root,"results/tables/disease_control_gse254300_cs/GSE254300_CS_paired_edgeR_full.tsv.gz")
 int_path<-file.path(root,"results/tables/disease_control_gse254300_cs/GSE254300_AIS_by_side_interaction_edgeR_full.tsv.gz")
 if(binding$status!="PASS"||binding$metadata_sha256!=specificity_hash(md_path)||
    binding$existing_A_counts_sha256!=specificity_hash(A_path)||
    binding$original_model_producer_sha256!=specificity_hash(producer_path)||
    binding$exporter_sha256!=specificity_hash(file.path(root,"scripts/export_specificity_counts.R"))||
    binding$original_CS_result_sha256!=specificity_hash(cs_path)||
    binding$original_interaction_result_sha256!=specificity_hash(int_path)||
    !isTRUE(binding$full_matrix_filterByExpr_verified)||!isTRUE(binding$raw_AIS_reused_count_values_exactly_equal))stop("Changed specificity input/source binding")
 seen<-character()
 for(f in binding$files) {
   if(length(f$path)!=1L||grepl("[/\\\\]",f$path)||f$path%in%c(".","..")||f$path%in%seen)stop("Invalid specificity input path")
   seen<-c(seen,f$path);p<-file.path(input,f$path)
   if(!file.exists(p)||nzchar(Sys.readlink(p))||file.info(p)$size!=f$bytes||specificity_hash(p)!=f$sha256)stop("Changed or missing specificity count input")
 }
 if(!setequal(seen,c("CS_joint_universe_counts.tsv.gz","AIS_joint_extra_counts.tsv.gz")))stop("Unbound specificity input")
 audit<-specificity_read(md_path);info<-audit[audit$cohort%in%c("AIS","congenital_scoliosis"),]
 stopifnot(nrow(info)==20L,all(table(info$cohort)==10L),all(table(info$cohort,info$side)==5L),all(table(info$patient_id)==2L))
 info$cohort_label<-ifelse(info$cohort=="AIS","AIS","CS");info$pair_id<-paste(info$cohort_label,info$patient_id,sep="_")
 info$side<-factor(info$side,levels=c("convex","concave"))
 expected_cs<-specificity_read(cs_path);expected_int<-specificity_read(int_path)
 stopifnot(nrow(expected_cs)==14450L,nrow(expected_int)==14967L,all(expected_cs$gene%in%expected_int$gene))
 A<-specificity_matrix(A_path);extra<-specificity_matrix(file.path(input,"AIS_joint_extra_counts.tsv.gz"));CS<-specificity_matrix(file.path(input,"CS_joint_universe_counts.tsv.gz"))
 counts_all<-specificity_reconstruct(A,extra,CS,expected_int$gene,info)
 set.seed(254305)
 # Original script61 joint tested universe and native TMM normalization.
 y_all<-calcNormFactors(DGEList(counts_all))
 original_qc<-specificity_read(file.path(root,"results/tables/disease_control_gse254300_cs/GSE254300_AIS_CS_joint_QC.tsv"))
 stopifnot(identical(original_qc$sample_id,colnames(counts_all)))
 specificity_close(y_all$samples$norm.factors,original_qc$normalization_factor,"joint native TMM factors",tolerance=2e-10)
 cs_idx<-info$cohort_label=="CS";info_cs<-droplevels(info[cs_idx,])
 y_cs<-calcNormFactors(DGEList(counts_all[expected_cs$gene,cs_idx,drop=FALSE]))
 design_cs<-model.matrix(~factor(patient_id)+side,info_cs);stopifnot(qr(design_cs)$rank==ncol(design_cs))
 y_cs<-estimateDisp(y_cs,design_cs,robust=TRUE);fit_cs<-glmQLFit(y_cs,design_cs,robust=TRUE)
 actual_cs<-topTags(glmQLFTest(fit_cs,coef="sideconcave"),n=Inf,sort.by="none")$table
 actual_cs$gene<-rownames(actual_cs);actual_cs$contrast<-"CS_concave_minus_convex";actual_cs$n_paired_patients<-5L
 info$side_concave<-as.integer(info$side=="concave");info$cohort_AIS<-as.integer(info$cohort_label=="AIS")
 design_int<-model.matrix(~factor(pair_id)+side_concave+side_concave:cohort_AIS,info)
 stopifnot(qr(design_int)$rank==ncol(design_int))
 y_int<-estimateDisp(y_all,design_int,robust=TRUE);fit_int<-glmQLFit(y_int,design_int,robust=TRUE)
 coefficient<-grep("side_concave:cohort_AIS",colnames(design_int),fixed=TRUE);stopifnot(length(coefficient)==1L)
 actual_int<-topTags(glmQLFTest(fit_int,coef=coefficient),n=Inf,sort.by="none")$table
 actual_int$gene<-rownames(actual_int);actual_int$contrast<-"AIS_minus_CS_difference_in_concave_minus_convex_effect"
 actual_int$n_AIS_pairs<-5L;actual_int$n_CS_pairs<-5L
 comparisons<-list(CS=specificity_compare(actual_cs,expected_cs,"CS"),interaction=specificity_compare(actual_int,expected_int,"interaction"))
 dir.create(run,recursive=TRUE,showWarnings=FALSE)
 write_gz<-function(x,name){con<-gzfile(file.path(run,name),"wt");on.exit(close(con));write.table(x,con,sep="\t",quote=FALSE,row.names=FALSE)}
 write_gz(actual_cs,"GSE254300_CS_paired_edgeR_full.tsv.gz");write_gz(actual_int,"GSE254300_AIS_by_side_interaction_edgeR_full.tsv.gz")
 normalization<-data.frame(model=c(rep("CS",10),rep("joint_interaction",20)),sample_id=c(colnames(y_cs),colnames(y_all)),normalization_factor=c(y_cs$samples$norm.factors,y_all$samples$norm.factors))
 write.table(normalization,file.path(run,"refitted_normalization_factors.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 report<-list(status="PASS",scope="Actual refits of both original script61 gene-level robust edgeR QL models from exact retained integer counts; no new filtering, model, disease-specificity inference or reference replacement",
 CS_genes_refitted=nrow(actual_cs),interaction_genes_refitted=nrow(actual_int),n_AIS_pairs=5L,n_CS_pairs=5L,
 CS_design="~ factor(patient_id) + side",interaction_design="~ factor(pair_id) + side_concave + side_concave:cohort_AIS",robust_dispersion=TRUE,robust_QL_fit=TRUE,
 normalization="Original model-specific native TMM; joint factors independently match original QC at unchanged arithmetic tolerance 2e-10",
 comparison_tolerance="Existing model-refit criterion 1e-7 * max(1,abs(reference)); exact genes/schema/counts, logFC signs and P/FDR<0.05 decisions",
 comparison=comparisons,CS_FDR05=sum(actual_cs$FDR<.05),interaction_FDR05=sum(actual_int$FDR<.05),
 refitter_sha256=specificity_hash(script),export_binding_sha256=specificity_hash(binding_path),existing_A_counts_sha256=specificity_hash(A_path),
 metadata_sha256=specificity_hash(md_path),original_model_producer_sha256=specificity_hash(producer_path),
 R_version=as.character(getRversion()),edgeR_version=as.character(packageVersion("edgeR")))
 jsonlite::write_json(report,file.path(run,"independent_specificity_refit.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
 message("PASS: original CS14450 and interaction14967 robust QL gene models refitted and independently compared")
 cat(normalizePath(run,mustWork=TRUE),"\n",sep="")
}
if(sys.nframe()==0L)specificity_main()
