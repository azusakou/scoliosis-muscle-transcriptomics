#!/usr/bin/env Rscript
# Contract: show the narrower genome-wide observation after the unchanged RNA
# gates, alongside the two fixed retrospective candidates and their uncertainty.
# R workflow, quantitative grid, 183 x 142 mm. No signature retraining.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(scales)})
source("scripts/lib_common.R")
rev <- "analysis"
rd <- function(n) read.delim(if(grepl("gz$",n))gzfile(file.path(rev,n))else file.path(rev,n),check.names=FALSE)
summary <- rd("bulk_gate_genomewide_exact_summary.tsv")
norm <- "locked_original_TMM"
s <- summary[summary$normalization==norm,]
a <- read.delim(gzfile(file.path(rev,"bulk_gate_inputs/original_A_gene_results.tsv.gz")),check.names=FALSE)
b <- rd(paste0("bulk_gate_observed_genes_",norm,".tsv.gz"))
genes <- read.delim(file.path(rev,"bulk_gate_inputs/common_gene_universe.tsv"))$gene
d <- data.frame(gene=genes,logFC_A=a$logFC[match(genes,a$gene)],logFC_B=b$logFC[match(genes,b$gene)])
stopifnot(nrow(d)==12995,!anyNA(d),nrow(s)==2)
null <- rd(paste0("bulk_gate_exact_null_",norm,".tsv"))
candidates <- rd("bulk_gate_ANGPT1_LEP_fixed_candidate_results.tsv")
candidates <- candidates[candidates$normalization==norm,]
patients <- rd("bulk_gate_fixed_candidate_patient_effects.tsv")
patients <- patients[patients$normalization==norm,]
write_tsv(d,file.path(rev,"tables/figure_S10_gate_gene_effects.tsv"))
write_tsv(null,file.path(rev,"tables/figure_S10_gate_exact_null.tsv"))
write_tsv(patients,file.path(rev,"tables/figure_S10_candidate_patient_effects.tsv"))
th <- theme_classic(base_size=7,base_family="Helvetica")+theme(
 plot.title=element_text(size=8,face="bold"),plot.subtitle=element_text(size=6.8),
 plot.tag=element_text(size=9,face="bold"),legend.title=element_text(size=6.5),
 legend.text=element_text(size=6),plot.margin=margin(5,7,5,5))
rho <- s$observed[s$endpoint=="spearman_logFC"]; agree <- s$observed[s$endpoint=="sign_concordance"]
p1 <- ggplot(d,aes(logFC_A,logFC_B))+geom_bin_2d(bins=50)+
 geom_hline(yintercept=0,colour="#A5ABB0",linewidth=.25)+geom_vline(xintercept=0,colour="#A5ABB0",linewidth=.25)+
 scale_fill_gradient(low="#DCE7EC",high="#315B77",trans="log10",name="Genes per bin")+
 labs(tag="a",title=NULL,subtitle=NULL,
  x="GSE254300 fitted log2 fold change",y="RNA-gate-passing GSE311180\nfitted log2 fold change")+th
p2 <- ggplot(null,aes(spearman_rho))+geom_histogram(bins=32,fill="#A8BEC9",colour="white",linewidth=.2)+
 geom_vline(xintercept=c(-rho,rho),colour="#B86545",linewidth=.6)+
 labs(tag="b",title=NULL,subtitle=NULL,
  x="Spearman correlation",y="Joint side configurations")+th
p3 <- ggplot(null,aes(sign_concordance))+geom_histogram(bins=32,fill="#A8BEC9",colour="white",linewidth=.2)+
 geom_vline(xintercept=c(1-agree,agree),colour="#B86545",linewidth=.6)+
 scale_x_continuous(labels=percent_format(accuracy=1))+
 labs(tag="c",title=NULL,subtitle=NULL,
  x="Common-gene sign agreement",y="Joint side configurations")+th
patients$gene <- factor(patients$gene,levels=c("ANGPT1","LEP"))
candidates$gene <- factor(candidates$gene,levels=c("ANGPT1","LEP"))
patients$display_x <- as.numeric(patients$gene)+rep(seq(-.15,.15,length.out=5),2)
p4 <- ggplot(patients,aes(display_x,paired_logCPM_difference))+
 geom_hline(yintercept=0,colour="#A5ABB0",linewidth=.3)+
 geom_point(size=1.7,colour="#315B77")+
 geom_errorbar(data=candidates,aes(x=as.numeric(gene),ymin=patient_t_CI95_lower,ymax=patient_t_CI95_upper),
  inherit.aes=FALSE,width=.28,linewidth=.6,colour="#B86545")+
 geom_point(data=candidates,aes(x=as.numeric(gene),y=patient_paired_logCPM_mean),inherit.aes=FALSE,
  shape=18,size=2.8,colour="#B86545")+
 scale_x_continuous(breaks=1:2,labels=c("ANGPT1","LEP"),limits=c(.6,2.4))+
 labs(tag="d",title=NULL,subtitle=NULL,
  x=NULL,y="Paired log2 CPM difference\n(concave minus convex)")+th
p <- (p1+p2)/(p3+p4)+plot_layout(guides="collect") & theme(legend.position="bottom")
dir.create(file.path(rev,"figures"),recursive=TRUE,showWarnings=FALSE)
save_pub_r(p,file.path(rev,"figures/figure_S10_gate_concordance"),183,142,600)
ragg::agg_png(file.path(rev,"figures/figure_S10_gate_concordance.png"),width=183,height=142,units="mm",res=180)
print(p);dev.off()
