#!/usr/bin/env Rscript
# Refit the original script72 diagnostic models from genuine retained counts.
# No new model, biological inference, target selection or reference replacement.

mixed_hash <- function(p) digest::digest(file=p,algo="sha256")
mixed_read <- function(p) read.delim(if(grepl("[.]gz$",p))gzfile(p) else p,check.names=FALSE,stringsAsFactors=FALSE)
mixed_key <- function(x, fields) do.call(paste,c(x[fields],sep="\034"))
mixed_close <- function(actual,expected,label,tolerance=1e-7) {
  if(length(actual)!=length(expected) || !identical(is.na(actual),is.na(expected)) ||
     any(!is.finite(actual)) || any(!is.finite(expected)) ||
     any(abs(actual-expected)>tolerance*pmax(1,abs(expected)))) stop("Mismatch: ",label)
}
mixed_check_files <- function(directory, files) {
  seen <- character()
  for(f in files) {
    name <- f$path
    if(length(name)!=1L || !nzchar(name) || grepl("[/\\\\]",name) || name %in% c(".","..") || name %in% seen)
      stop("Unsafe or duplicate mixed-model input artifact")
    seen <- c(seen,name); p <- file.path(directory,name)
    if(!file.exists(p) || nzchar(Sys.readlink(p)) || file.info(p)$size!=f$bytes || mixed_hash(p)!=f$sha256)
      stop("Changed or missing mixed-model input: ",name)
  }
  invisible(seen)
}
mixed_validate_metadata <- function(metadata) {
  required <- c("cell_id","cell_type","sample_id","patient_id","side","nCount_RNA","nCount_SoupX","subsets_Mito_percent")
  if(!all(required %in% names(metadata)) || anyNA(metadata[required]) ||
     anyDuplicated(metadata$cell_id) || any(!metadata$side %in% c("convex","concave")) ||
     any(!is.finite(metadata$nCount_RNA)) || any(!is.finite(metadata$nCount_SoupX)) ||
     any(metadata$nCount_RNA<=0) || any(metadata$nCount_SoupX<=0) ||
     any(!is.finite(metadata$subsets_Mito_percent))) stop("Invalid mixed-model cell metadata")
  invisible(TRUE)
}
mixed_sidecar_values <- function(pair,metadata,counts,availability) {
  a <- availability[availability$assay=="SoupX" & availability$cell_type==pair$cell_type & availability$gene==pair$gene,,drop=FALSE]
  cells <- metadata$cell_id[metadata$cell_type==pair$cell_type]
  if(nrow(a)!=1L || !isTRUE(a$available) || a$n_cells!=length(cells) ||
     a$all_gene_total_source_field!="nCount_SoupX" || a$zero_encoding!="every_model_cell_explicit_including_zero")
    stop("Missing or invalid model-specific availability: ",pair$gene)
  z <- counts[counts$assay=="SoupX" & counts$cell_type==pair$cell_type & counts$gene==pair$gene,,drop=FALSE]
  if(nrow(z)!=length(cells) || anyDuplicated(z$cell_id) || !setequal(z$cell_id,cells) ||
     any(!is.finite(z$UMI)) || any(z$UMI<0) || any(z$UMI!=round(z$UMI)) || sum(z$UMI>0)!=a$n_nonzero)
    stop("Incomplete or invalid explicit sidecar cell coverage: ",pair$gene)
  z$UMI[match(cells,z$cell_id)]
}

mixed_refit <- function(root,input,extra,run,script) {
  suppressPackageStartupMessages(library(lme4))
  root <- normalizePath(root,mustWork=TRUE); input <- normalizePath(input,mustWork=TRUE)
  extra <- normalizePath(extra,mustWork=TRUE)
  if(file.exists(run)) stop("Mixed-model output already exists; use a new directory")
  binding_path <- file.path(root,"metadata/audit/retained_cell_counts_export.json")
  metadata_path <- file.path(root,"results/tables/annotation/annotated_cell_metadata.tsv.gz")
  de_path <- file.path(root,"results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz")
  producer_path <- file.path(root,"scripts/72_model_sensitivities.R")
  fixed_binding <- jsonlite::fromJSON(binding_path,simplifyVector=FALSE)
  if(fixed_binding$status!="PASS" || mixed_hash(metadata_path)!=fixed_binding$metadata_sha256 ||
     mixed_hash(file.path(root,"scripts/export_retained_cell_counts.R"))!=fixed_binding$exporter_sha256)
    stop("Changed fixed sparse-count source binding")
  mixed_check_files(input,fixed_binding$files)
  extra_binding_path <- file.path(extra,"mixed_model_extra_export.json")
  if(!file.exists(extra_binding_path)) stop("Missing model sidecar binding; unexported genes cannot be zero-filled")
  extra_binding <- jsonlite::fromJSON(extra_binding_path,simplifyVector=FALSE)
  if(extra_binding$status!="PASS" || extra_binding$metadata_sha256!=mixed_hash(metadata_path) ||
     extra_binding$fixed_target_binding_sha256!=mixed_hash(binding_path) ||
     extra_binding$input_object_sha256!=fixed_binding$input_object_sha256 ||
     extra_binding$selected_original_pseudobulk_sha256!=mixed_hash(de_path) ||
     extra_binding$original_model_producer_sha256!=mixed_hash(producer_path) ||
     extra_binding$exporter_sha256!=mixed_hash(file.path(root,"scripts/export_cell_mixed_sidecar.R")))
    stop("Changed mixed-model sidecar/source binding")
  mixed_check_files(extra,extra_binding$files)
  metadata <- mixed_read(metadata_path); mixed_validate_metadata(metadata)
  stopifnot(nrow(metadata)==14846L,length(unique(metadata$patient_id))==2L,
            all(metadata$passes_qc),all(metadata$scDblFinder.class=="singlet"))
  de <- mixed_read(de_path)
  selected <- do.call(rbind,lapply(split(de,de$cell_type),function(z) head(z[order(abs(z$mean_donor_effect),decreasing=TRUE),],10)))
  stopifnot(nrow(selected)==50L,!anyDuplicated(selected[,c("cell_type","gene")]))
  target_availability <- mixed_read(file.path(input,"cell_target_availability.tsv"))
  covered <- target_availability$gene[target_availability$assay=="SoupX" & target_availability$available]
  extra_availability <- mixed_read(file.path(extra,"mixed_model_extra_availability.tsv"))
  extra_counts <- mixed_read(file.path(extra,"mixed_model_extra_counts.tsv.gz"))
  missing <- selected[!selected$gene %in% covered,c("cell_type","gene"),drop=FALSE]
  if(!setequal(mixed_key(missing,c("cell_type","gene")),mixed_key(extra_availability,c("cell_type","gene"))) ||
     anyDuplicated(extra_availability[,c("assay","cell_type","gene")])) stop("Sidecar target set changed")
  # Zero-fill is allowed solely after a true assay/gene availability assertion.
  counts <- mixed_read(file.path(input,"cell_target_counts.tsv.gz"))
  counts <- counts[counts$assay=="SoupX" & counts$gene %in% selected$gene,,drop=FALSE]
  stopifnot(!anyDuplicated(counts[,c("cell_id","gene")]),all(counts$cell_id %in% metadata$cell_id),
            all(is.finite(counts$UMI)),all(counts$UMI>0),all(counts$UMI==round(counts$UMI)))
  for(i in seq_len(nrow(missing))) mixed_sidecar_values(missing[i,],metadata,extra_counts,extra_availability)
  lock <- mixed_read(file.path(root,"environment.lock.tsv"))
  lock_lme4 <- as.character(lock[lock[[1]]=="lme4",2])
  if(length(lock_lme4)!=1L || as.character(packageVersion("lme4"))!=lock_lme4)
    stop("lme4 runtime differs from the published model runtime")
  set.seed(720014)
  models <- list(); details <- list()
  for(i in seq_len(nrow(selected))) {
    z <- selected[i,]; cells <- which(metadata$cell_type==z$cell_type); ids <- metadata$cell_id[cells]
    if(z$gene %in% covered) {
      umi <- numeric(length(cells)); source <- counts[counts$gene==z$gene & counts$cell_id %in% ids,,drop=FALSE]
      umi[match(source$cell_id,ids)] <- source$UMI
      source_role <- "fixed_165_target_sparse_counts_explicit_available_gene"
    } else {
      umi <- mixed_sidecar_values(z,metadata,extra_counts,extra_availability)
      source_role <- "same_object_missing_gene_dense_lineage_sidecar"
    }
    stopifnot(all(umi<=metadata$nCount_SoupX[cells]))
    dat <- data.frame(value=log1p(umi/metadata$nCount_SoupX[cells]*1e4),
      side=factor(metadata$side[cells],levels=c("convex","concave")),
      patient=metadata$patient_id[cells],log_umi=log10(metadata$nCount_RNA[cells]+1),
      mt=metadata$subsets_Mito_percent[cells])
    # Verbatim formula, ML/REML choice and controls from original script72.
    fit <- suppressWarnings(lmer(value~side+log_umi+mt+(1|patient),dat,REML=FALSE,
      control=lmerControl(check.nlev.gtr.1="ignore",check.nobs.vs.nlev="ignore",check.nobs.vs.rankZ="ignore")))
    co <- coef(summary(fit))["sideconcave",]
    models[[i]] <- data.frame(cell_type=z$cell_type,gene=z$gene,pseudobulk_effect=z$mean_donor_effect,
      cell_mixed_effect=co["Estimate"],cell_mixed_se=co["Std. Error"],
      sign_concordant=sign(z$mean_donor_effect)==sign(co["Estimate"]),singular_fit=isSingular(fit),
      inferential_status="diagnostic_only_two_patient_random_effect")
    details[[i]] <- data.frame(cell_type=z$cell_type,gene=z$gene,n_cells=nrow(dat),
      n_patients=length(unique(dat$patient)),n_nonzero=sum(umi>0),input_role=source_role,
      convergence_messages=paste(fit@optinfo$conv$lme4$messages,collapse=";"))
  }
  actual <- do.call(rbind,models)
  expected_path <- file.path(root,"results/tables/robustness/pseudobulk_vs_cell_mixed_model.tsv")
  expected <- mixed_read(expected_path)
  stopifnot(identical(mixed_key(actual,c("cell_type","gene")),mixed_key(expected,c("cell_type","gene"))))
  for(field in c("pseudobulk_effect","cell_mixed_effect","cell_mixed_se"))
    mixed_close(actual[[field]],expected[[field]],field)
  for(field in c("sign_concordant","singular_fit","inferential_status"))
    if(!identical(actual[[field]],expected[[field]])) stop("Exact model flag/status changed: ",field)
  summary <- aggregate(sign_concordant~cell_type,actual,mean)
  expected_summary <- mixed_read(file.path(root,"results/tables/robustness/pseudobulk_vs_cell_mixed_model_summary.tsv"))
  stopifnot(identical(summary$cell_type,expected_summary$cell_type),
            identical(summary$sign_concordant,expected_summary$sign_concordant))
  dir.create(run,recursive=TRUE,showWarnings=FALSE)
  write.table(actual,file.path(run,"pseudobulk_vs_cell_mixed_model.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(summary,file.path(run,"pseudobulk_vs_cell_mixed_model_summary.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  write.table(do.call(rbind,details),file.path(run,"model_input_and_convergence.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
  report <- list(status="PASS",scope="Actual ML refits of original diagnostic cell mixed models from authentic target counts and full-assay QC metadata; no new patient-level inference, discovery or upstream cell reconstruction",
    models_refitted=nrow(actual),sign_concordant=sum(actual$sign_concordant),singular_models=sum(actual$singular_fit),
    formula="value ~ side + log_umi + mt + (1 | patient)",REML=FALSE,
    control=list(check.nlev.gtr.1="ignore",check.nobs.vs.nlev="ignore",check.nobs.vs.rankZ="ignore"),
    transformation="SoupX log1p(UMI / full_assay_UMI * 10000); log_umi = log10(raw_RNA_all_gene_UMI + 1); mt = raw mitochondrial percent",
    comparison_tolerance="Existing model-refit criterion: 1e-7 * max(1,abs(reference)); exact cell/gene IDs, signs, singularity and status",
    max_absolute_difference=lapply(c("pseudobulk_effect","cell_mixed_effect","cell_mixed_se"),function(f)list(field=f,max_abs=max(abs(actual[[f]]-expected[[f]])))),
    refitter_sha256=mixed_hash(script),original_model_producer_sha256=mixed_hash(producer_path),
    fixed_target_binding_sha256=mixed_hash(binding_path),extra_binding_sha256=mixed_hash(extra_binding_path),
    metadata_sha256=mixed_hash(metadata_path),original_pseudobulk_sha256=mixed_hash(de_path),reference_model_sha256=mixed_hash(expected_path),
    R_version=as.character(getRversion()),lme4_version=as.character(packageVersion("lme4")),Matrix_version=as.character(packageVersion("Matrix")))
  jsonlite::write_json(report,file.path(run,"independent_cell_mixed_refit.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
  run
}

mixed_main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  allowed <- c("--root","--input-dir","--extra-dir","--output-dir")
  flags <- sub("=.*$","",args)
  if(any(!flags %in% allowed) || any(!grepl("=.+$",args)) || anyDuplicated(flags)) stop("Unknown, empty or duplicate CLI flag")
  value <- function(flag,default) {z<-args[startsWith(args,paste0(flag,"="))];if(length(z))sub(paste0("^",flag,"="),"",z) else default}
  script <- normalizePath(sub("^--file=","",commandArgs(FALSE)[grep("^--file=",commandArgs(FALSE))][1]),mustWork=TRUE)
  root <- value("--root",dirname(dirname(script)))
  input <- value("--input-dir",file.path(root,"data/retained"))
  extra <- value("--extra-dir",file.path(root,"data/retained/mixed_model_extra"))
  logs <- file.path(root,"logs/cell_mixed_model_replay");dir.create(logs,recursive=TRUE,showWarnings=FALSE)
  run <- value("--output-dir",tempfile("run-",tmpdir=logs))
  result <- mixed_refit(root,input,extra,run,script)
  message("PASS: all original mixed models refitted and exact sign/singularity flags matched the immutable references")
  cat(normalizePath(result,mustWork=TRUE),"\n",sep="")
}
if(sys.nframe()==0L) mixed_main()
