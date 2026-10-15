#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(lme4); library(harmony) })
source("scripts/lib_common.R")
set.seed(720014)
options(future.globals.maxSize = 6 * 1024^3)
future::plan("sequential")

obj <- readRDS("data/processed/scrna/scrna_annotated.rds")
obj <- JoinLayers(obj, assay="SoupX")
expr <- LayerData(obj, assay="SoupX", layer="data")
de <- read.delim(gzfile("results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz"),check.names=FALSE)
selected <- do.call(rbind,lapply(split(de,de$cell_type),function(z) head(z[order(abs(z$mean_donor_effect),decreasing=TRUE),],10)))
mixed <- list()
for(i in seq_len(nrow(selected))){
  z <- selected[i,]; cells <- rownames(obj@meta.data)[obj$cell_type==z$cell_type]
  if(!z$gene %in% rownames(expr)) next
  dat <- data.frame(value=as.numeric(expr[z$gene,cells]),side=factor(obj$side[match(cells,colnames(obj))],levels=c("convex","concave")),
                    patient=obj$patient_id[match(cells,colnames(obj))],log_umi=log10(obj$nCount_RNA[match(cells,colnames(obj))]+1),
                    mt=obj$subsets_Mito_percent[match(cells,colnames(obj))])
  fit <- suppressWarnings(lmer(value~side+log_umi+mt+(1|patient),dat,REML=FALSE,
                               control=lmerControl(check.nlev.gtr.1="ignore",check.nobs.vs.nlev="ignore",check.nobs.vs.rankZ="ignore")))
  co <- coef(summary(fit))["sideconcave",]
  mixed[[length(mixed)+1]] <- data.frame(cell_type=z$cell_type,gene=z$gene,pseudobulk_effect=z$mean_donor_effect,
    cell_mixed_effect=co["Estimate"],cell_mixed_se=co["Std. Error"],sign_concordant=sign(z$mean_donor_effect)==sign(co["Estimate"]),
    singular_fit=isSingular(fit),inferential_status="diagnostic_only_two_patient_random_effect")
}
mixed <- do.call(rbind,mixed)
write_tsv(mixed,"results/tables/robustness/pseudobulk_vs_cell_mixed_model.tsv")
write_tsv(aggregate(sign_concordant~cell_type,mixed,mean),"results/tables/robustness/pseudobulk_vs_cell_mixed_model_summary.tsv")

# Quantify whether batch-corrected embeddings retain or erase side structure.
# These embeddings never enter DE or program scoring.
obj <- RunHarmony(obj,group.by.vars="sample_id",reduction.use="pca",reduction.save="harmony",plot_convergence=FALSE,verbose=FALSE)
embedding_metrics <- function(mat,label,meta){
  pcs <- seq_len(min(20,ncol(mat))); out <- list()
  for(j in pcs){
    d <- data.frame(y=mat[,j],patient=meta$patient_id,side=meta$side,sample=meta$sample_id)
    side_fit <- lm(y~factor(patient)+factor(side),d); sample_fit <- lm(y~factor(sample),d)
    # Sequential side SS after patient adjustment divided by total SS is
    # semi-partial R2, not partial R2 (which uses side SS + residual SS).
    out[[j]] <- data.frame(method=label,component=j,side_semipartial_r2=anova(side_fit)["factor(side)","Sum Sq"]/sum(anova(side_fit)$"Sum Sq"),
                           sample_r2=summary(sample_fit)$r.squared)
  }
  do.call(rbind,out)
}
metrics <- rbind(embedding_metrics(Embeddings(obj,"pca"),"unintegrated_PCA",obj@meta.data),
                 embedding_metrics(Embeddings(obj,"harmony"),"Harmony",obj@meta.data))

# Seurat reciprocal-PCA sensitivity.
rpca <- obj
rpca[["SoupX"]] <- split(rpca[["SoupX"]],f=rpca$sample_id)
rpca <- NormalizeData(rpca,verbose=FALSE); rpca <- FindVariableFeatures(rpca,nfeatures=2000,verbose=FALSE)
rpca <- ScaleData(rpca,verbose=FALSE); rpca <- RunPCA(rpca,npcs=30,verbose=FALSE)
rpca <- IntegrateLayers(rpca,method=RPCAIntegration,orig.reduction="pca",new.reduction="integrated.rpca",verbose=FALSE)
metrics <- rbind(metrics,embedding_metrics(Embeddings(rpca,"integrated.rpca"),"Seurat_RPCA",rpca@meta.data))
write_tsv(metrics,"results/tables/robustness/integration_embedding_sensitivity.tsv")
write_tsv(aggregate(cbind(side_semipartial_r2,sample_r2)~method,metrics,mean),"results/tables/robustness/integration_embedding_sensitivity_summary.tsv")
message("Mixed-model and integration sensitivities complete.")
