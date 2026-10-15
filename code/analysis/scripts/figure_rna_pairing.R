#!/usr/bin/env Rscript
# Complete RNA profile comparisons, metadata-pair outlines and actual coverage.
# Descriptive only. No anatomical/clinical identity or biological replication claim.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(scales)})
source("scripts/lib_common.R")
rev <- "analysis"
pairs <- read.delim(file.path(rev, "audit_gse311180_rna_all_pairs.tsv"), check.names=FALSE)
expected <- read.delim(analysis_input(file.path(rev, "audit_gse311180_rna_expected_pairs.tsv")), check.names=FALSE)
meta <- read.delim("metadata/audit/gse311180_requant_manifest.tsv", check.names=FALSE)
stopifnot(nrow(meta)==14, !anyDuplicated(meta$sample_id), nrow(pairs)==91,
          nrow(expected)==7, all(is.finite(pairs$relatedness)),
          all(pairs$n>=0), all(pairs$ibs0>=0))
meta <- meta[order(meta$patient_id, factor(meta$side,levels=c("concave","convex"))),]
samples <- meta$sample_id
labels <- setNames(paste0(meta$sample_id," (",ifelse(meta$side=="concave","Cc","Cv"),")"),samples)
forward <- data.frame(row=pairs$sample_a,col=pairs$sample_b,relatedness=pairs$relatedness,
                      n=pairs$n,ibs0=pairs$ibs0)
reverse <- forward;reverse$row <- forward$col;reverse$col <- forward$row
grid <- rbind(forward,reverse,data.frame(row=samples,col=samples,relatedness=NA_real_,n=NA,ibs0=NA))
grid$row <- factor(grid$row,levels=rev(samples));grid$col <- factor(grid$col,levels=samples)
grid$expected_pair <- as.character(grid$row)!=as.character(grid$col) &
  meta$patient_id[match(grid$row,samples)]==meta$patient_id[match(grid$col,samples)]
stopifnot(nrow(grid)==196, sum(grid$expected_pair)==14,
          !anyDuplicated(paste(grid$row,grid$col)))
range <- range(c(-.1,1,grid$relatedness),na.rm=TRUE)
th <- theme_classic(base_size=7,base_family="Helvetica") + theme(
  plot.tag=element_text(face="bold",size=9),plot.title=element_text(face="bold",size=8),
  plot.subtitle=element_text(size=6.5),legend.title=element_text(size=6.5),
  legend.text=element_text(size=6),plot.margin=margin(5,5,5,5))
a <- ggplot(grid,aes(col,row,fill=relatedness)) + geom_tile(colour="white",linewidth=.2) +
  geom_tile(data=grid[grid$expected_pair,],fill=NA,colour="#222222",linewidth=.55) +
  geom_point(data=grid[!is.na(grid$n)&grid$n<1000,],shape=4,size=1.1,colour="#C56A48") +
  scale_fill_gradient(low="#F7F7F2",high="#355C8A",limits=range,na.value="#DDDDDD",
                      name="RNA profile\nrelatedness") +
  scale_x_discrete(labels=labels,drop=FALSE) + scale_y_discrete(labels=labels,drop=FALSE) +
  coord_equal() + labs(tag="a",x=NULL,y=NULL,title="All 91 library comparisons",
    subtitle="Black outlines: seven expected metadata pairs\nCrosses: fewer than 1,000 shared callable sites") + th +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6),
        axis.text.y=element_text(size=6),axis.ticks=element_blank(),
        axis.line=element_blank(),legend.position="bottom")
expected$patient_id <- factor(expected$patient_id,levels=rev(unique(meta$patient_id)))
expected$coverage_label <- sprintf("n = %s; IBS0 = %s",
                                  format(expected$shared_callable_sites,big.mark=",",trim=TRUE),
                                  format(expected$ibs0,trim=TRUE))
expected$verdict <- ifelse(expected$interpretation=="consistent_RNA_profiles","Meets RNA gates",
                           ifelse(expected$interpretation=="inconclusive_coverage","Low coverage","Does not meet gates"))
b <- ggplot(expected,aes(relatedness,patient_id,colour=verdict)) +
  geom_vline(xintercept=.9,linetype="dashed",linewidth=.3,colour="#888888") +
  geom_point(size=2.1) +
  geom_text(aes(x=.85,label=coverage_label),hjust=0,vjust=1.8,size=2.2,colour="#333333") +
  scale_colour_manual(values=c("Meets RNA gates"="#355C8A","Low coverage"="#888888",
                               "Does not meet gates"="#C56A48"),name=NULL) +
  scale_x_continuous(limits=c(.85,1.005),breaks=c(.85,.9,.95,1),expand=expansion(mult=c(.025,.05))) +
  labs(tag="b",x="RNA profile relatedness",y=NULL,title="Expected-pair results",
       subtitle="Predeclared consistency gates; no confidence intervals\nPairing check does not verify anatomical side") + th +
  theme(legend.position="bottom",legend.text=element_text(size=6))
figure <- a+b+plot_layout(widths=c(1.6,1))
dir.create(file.path(rev,"figures"),recursive=TRUE,showWarnings=FALSE)
save_pub_r(figure,file.path(rev,"figures/figure_S8_RNA_pairing"),183,150,600)
ragg::agg_png(file.path(rev,"figures/figure_S8_RNA_pairing.png"),width=183,height=150,units="mm",res=180)
print(figure);dev.off()
write_tsv(grid,file.path(rev,"tables/figure_S8_RNA_pairing_all_comparisons.tsv"))
write_tsv(expected,file.path(rev,"tables/figure_S8_RNA_expected_pairs.tsv"))
