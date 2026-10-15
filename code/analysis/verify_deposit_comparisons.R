#!/usr/bin/env Rscript
# R's built-in tie-ranked cor is independent of the Python producer's rank code.
raw<-read.delim(gzfile("analysis/bulk_followup_full_gene_TPM_counts.tsv.gz"),check.names=FALSE)
dep<-read.delim(gzfile("analysis/followup_inputs/GSE311180_deposited_TPM.txt.gz"),check.names=FALSE)
dep<-dep[!is.na(dep$gene_name) & dep$gene_name!="",]
genes<-sort(intersect(raw$gene,dep$gene_name));stopifnot(length(genes)==22816L)
raw<-raw[match(genes,raw$gene),];dep<-dep[match(genes,dep$gene_name),]
tab<-read.delim("analysis/audit_gse311180_all_column_expression_comparisons.tsv",check.names=FALSE)
stopifnot(nrow(tab)==336L,all(tab$shared_named_genes==length(genes)))
for(i in seq_len(nrow(tab))) {
  x<-raw[[paste0(tab$raw_sample[i],"_TPM")]];y<-dep[[tab$deposited_column[i]]]
  stopifnot(abs(cor(x,y,method="spearman")-tab$spearman_TPM[i])<1e-12,
            abs(cor(log1p(x),log1p(y))-tab$pearson_log1p_TPM[i])<1e-12)
}
cat("PASS: all336 Python-produced expression comparisons independently recomputed in R; no labels authenticated or changed.\n")
