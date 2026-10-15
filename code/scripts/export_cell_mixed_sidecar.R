#!/usr/bin/env Rscript
# Export only original script72 model inputs missing from the fixed 165-target
# sparse-count release. No target, cell label, formula or reference fit changes.
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1L || length(args) > 2L) stop("Usage: Rscript scripts/export_cell_mixed_sidecar.R reconstructed_object.rds [NEW_OUTPUT_DIRECTORY]")
script <- normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]), mustWork=TRUE)
root <- dirname(dirname(script))
input <- normalizePath(args[1], mustWork=TRUE)
out <- if(length(args)>1L) args[2] else file.path(root,"data/retained/mixed_model_extra")
if(file.exists(out)) stop("Sidecar output already exists; use a new directory")
hash <- function(p) digest::digest(file=p,algo="sha256")
rd <- function(p) read.delim(if(grepl("[.]gz$",p))gzfile(p) else p, check.names=FALSE, stringsAsFactors=FALSE)
write_tsv <- function(x,p) {
  con <- if(grepl("[.]gz$",p)) gzfile(p,"wt") else file(p,"wt")
  on.exit(close(con)); write.table(x,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
metadata_path <- file.path(root,"results/tables/annotation/annotated_cell_metadata.tsv.gz")
binding_path <- file.path(root,"metadata/audit/retained_cell_counts_export.json")
de_path <- file.path(root,"results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz")
availability_path <- file.path(root,"data/retained/cell_target_availability.tsv")
binding <- jsonlite::fromJSON(binding_path,simplifyVector=FALSE)
stopifnot(binding$status=="PASS",hash(input)==binding$input_object_sha256,
          hash(metadata_path)==binding$metadata_sha256)
metadata <- rd(metadata_path); availability <- rd(availability_path); de <- rd(de_path)
selected <- do.call(rbind,lapply(split(de,de$cell_type),function(z) head(z[order(abs(z$mean_donor_effect),decreasing=TRUE),],10)))
stopifnot(nrow(selected)==50L,!anyDuplicated(selected[,c("cell_type","gene")]))
covered <- availability$gene[availability$assay=="SoupX" & availability$available]
missing <- selected[!selected$gene %in% covered,c("cell_type","gene"),drop=FALSE]
if(!nrow(missing)) stop("No missing model inputs require a sidecar")
object <- readRDS(input); md <- object@meta.data
stopifnot(nrow(md)==14846L, identical(rownames(md),metadata$cell_id),
          !anyDuplicated(metadata$cell_id))
for(f in c("sample_id","patient_id","side","cell_type","passes_qc","scDblFinder.class"))
  if(!identical(as.character(md[[f]]),as.character(metadata[[f]]))) stop("Metadata identity mismatch: ",f)
for(f in c("nCount_RNA","nCount_SoupX"))
  if(!identical(as.numeric(md[[f]]),as.numeric(metadata[[f]]))) stop("All-gene total mismatch: ",f)
stopifnot(all(abs(md$subsets_Mito_percent-metadata$subsets_Mito_percent)<=2e-12*pmax(1,abs(metadata$subsets_Mito_percent))))
object <- JoinLayers(object,assay="SoupX")
m <- LayerData(object,assay="SoupX",layer="counts")
normalized <- LayerData(object,assay="SoupX",layer="data")
stopifnot(identical(colnames(m),metadata$cell_id),
          all(Matrix::colSums(m)==metadata$nCount_SoupX),
          all(missing$gene %in% rownames(m)),all(missing$gene %in% rownames(normalized)))
counts <- list(); available <- list(); max_normalization_error <- 0
for(i in seq_len(nrow(missing))) {
  z <- missing[i,]; ix <- which(metadata$cell_type==z$cell_type); cells <- metadata$cell_id[ix]
  umi <- as.numeric(m[z$gene,cells])
  stopifnot(length(ix)>0L,all(is.finite(umi)),all(umi>=0),all(umi==round(umi)))
  actual <- log1p(umi/metadata$nCount_SoupX[ix]*1e4)
  stored <- as.numeric(normalized[z$gene,cells]); error <- max(abs(actual-stored))
  stopifnot(error <= 2e-10*max(1,abs(stored)))
  max_normalization_error <- max(max_normalization_error,error)
  # Dense within this original lineage: every zero is explicitly exported.
  counts[[i]] <- data.frame(assay="SoupX",cell_type=z$cell_type,gene=z$gene,
                           cell_id=cells,UMI=umi)
  available[[i]] <- data.frame(assay="SoupX",cell_type=z$cell_type,gene=z$gene,
    available=TRUE,n_cells=length(cells),n_nonzero=sum(umi>0),
    all_gene_total_source_field="nCount_SoupX",zero_encoding="every_model_cell_explicit_including_zero")
}
dir.create(out,recursive=TRUE,showWarnings=FALSE)
names <- c("mixed_model_extra_counts.tsv.gz","mixed_model_extra_availability.tsv")
write_tsv(do.call(rbind,counts),file.path(out,names[1]))
write_tsv(do.call(rbind,available),file.path(out,names[2]))
files <- lapply(names,function(n)list(path=n,bytes=unname(file.info(file.path(out,n))$size),sha256=hash(file.path(out,n))))
report <- list(status="PASS",scope="Actual missing original-script72 model counts; all lineage cells exported including real zero counts; no model inference or target change",
 input_object_role="same_independently_reconstructed_annotated_singlets_as_fixed_165_target_export",
 input_object_sha256=hash(input),exporter_sha256=hash(script),
 metadata_sha256=hash(metadata_path),fixed_target_binding_sha256=hash(binding_path),
 selected_original_pseudobulk_sha256=hash(de_path),original_model_producer_sha256=hash(file.path(root,"scripts/72_model_sensitivities.R")),
 n_missing_model_gene_lineage_pairs=nrow(missing),missing_model_pairs=missing,
 source_data_layer_normalization_max_abs_error=max_normalization_error,
 transform="Seurat SoupX LogNormalize: log1p(UMI / full_assay_UMI * 10000)",
 files=files,R_version=as.character(getRversion()),Seurat_version=as.character(packageVersion("Seurat")))
jsonlite::write_json(report,file.path(out,"mixed_model_extra_export.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
message("Missing model input export PASS: ",nrow(missing)," gene-lineage pair(s), ",sum(vapply(counts,nrow,integer(1)))," explicit cell records")
