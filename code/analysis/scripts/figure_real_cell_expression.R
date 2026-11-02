#!/usr/bin/env Rscript
# Contract: reveal sparse nominated signals and correction dependence in all
# four libraries and nine lineages; descriptive grid, no biological replication.
# R only. 183 x 190 mm; shared before/after scales; editable PDF/SVG and previews.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(scales)})
source("scripts/lib_common.R")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>1 || (length(args) && !startsWith(args,"--input=")))
 stop("Usage: Rscript analysis/scripts/figure_real_cell_expression.R [--input=summary.tsv]")
input <- if(length(args))substring(args,9) else "analysis/signature_expression_by_library_lineage_before_after.tsv"
stopifnot(nzchar(input),file.exists(input))
d <- read.delim(input,check.names=FALSE)
genes <- c("IGKC","ACTA1","MYL2","DES","CFD","APOD","POSTN","PRG4","MFAP5","ESR1")
short <- c("fibro-adipogenic progenitor"="FAP","endothelial cell"="Endothelial",
 "macrophage/monocyte"="Mac/Mono","muscle stem/progenitor cell"="MuSC",
 "T/NK cell"="T/NK","B/plasma cell"="B/plasma","pericyte"="Pericyte",
 "smooth-muscle-like cell"="SMC-like","Schwann/perineurial cell"="Schwann/perineurial")
d <- d[d$gene %in% genes,]
stopifnot(nrow(d)==720,!anyNA(short[d$cell_type]))
d$gene <- factor(d$gene,levels=genes)
d$lineage <- factor(unname(short[d$cell_type]),levels=rev(unname(short)))
libs <- c("H01045__concave","H01045__convex","H66206__concave","H66206__convex")
limit <- max(d$mean_log1p_cp10k,na.rm=TRUE)
write_tsv(d,"analysis/tables/figure_S7_real_cell_expression.tsv")
out <- "analysis/figures";dir.create(out,recursive=TRUE,showWarnings=FALSE)
for (assay in c("RNA","SoupX")) {
 z <- d[d$assay==assay,]
 plots <- lapply(seq_along(libs),function(i){
  s <- z[z$sample_id==libs[i],];n <- sum(s$n_cells[s$gene==genes[1]])
  ggplot(s,aes(gene,lineage)) +
   geom_tile(fill="white",colour="#E1E5E8",linewidth=.25) +
   geom_point(data=s[s$n_cells>0&s$n_detected>0,],
    aes(size=detection_fraction,colour=mean_log1p_cp10k)) +
   geom_point(data=s[s$n_cells==0,],shape=4,size=1.8,colour="#555555",stroke=.4) +
   scale_size_area(max_size=4.2,limits=c(0,1),breaks=c(.01,.1,.5,1),labels=percent_format(accuracy=1),
     name="Detected cells (%)") +
   scale_colour_gradientn(colours=c("#D4E4EC","#6892AE","#244B70"),limits=c(0,limit),
     name="Mean log1p CP10k",guide=guide_colourbar(barwidth=unit(25,"mm"),barheight=unit(2.5,"mm"))) +
   scale_x_discrete(drop=FALSE)+scale_y_discrete(drop=FALSE) +
   labs(tag=letters[i+if(assay=="RNA")0 else 4],x=NULL,y=NULL,
    title=NULL,subtitle=NULL,
    caption=paste0(sub("__", " | ", libs[i]), "; n = ", format(n,big.mark=","))) +
   theme_minimal(base_size=7,base_family="Helvetica") +
   theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1,colour="#252A31",face="italic"),
    axis.text.y=element_text(colour="#252A31"),plot.title=element_text(face="bold",size=8),
    plot.subtitle=element_text(size=7),plot.tag=element_text(face="bold",size=9),
    legend.title=element_text(size=7),legend.text=element_text(size=6.5),
    plot.margin=margin(6,4,6,4))
 })
 p <- wrap_plots(plots,ncol=2,guides="collect") + plot_layout(guides="collect") & theme(legend.position="bottom")
 stem <- paste0("figure_S7_",if(assay=="RNA")"before" else "after")
 save_pub_r(p,file.path(out,stem),183,190,600)
 ragg::agg_png(file.path(out,paste0(stem,".png")),width=183,height=190,units="mm",res=180)
 print(p);dev.off()
}
