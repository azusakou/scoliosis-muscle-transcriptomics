suppressPackageStartupMessages({library(Matrix); library(ggplot2)})
data_root <- Sys.getenv('AIS_DATA_ROOT','.')
out <- Sys.getenv('AIS_OUTPUT','computed')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
rd <- function(rel) read.delim(file.path(data_root,rel),check.names=FALSE,stringsAsFactors=FALSE)
wr <- function(x,n) write.table(x,file.path(out,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
cells <- rd('data/cell_counts/cells.tsv.gz')
nonzero <- rd('data/cell_counts/nonzero_counts.tsv.gz')
sig <- rd('results/tables/bulk/frozen_multicellular_signature.tsv')
ann <- read.delim('config/annotation_cluster_map.tsv',check.names=FALSE,stringsAsFactors=FALSE)
eff <- rd('results/tables/pseudobulk/scrna_donor_specific_effects.tsv.gz')
det <- rd('results/tables/bulk/frozen_multicellular_signature_projection_genes.tsv')
other <- rd('results/tables/bulk_gse311180_requant/GSE311180_label_blind_projection_logCPM.tsv.gz')$gene
corrected_full <- rd('analysis/length_aware/full_scaledTPM_projection_genes.tsv')$gene
corrected_gate <- rd('analysis/length_aware/gate_scaledTPM_projection_genes.tsv')$gene
markers <- c('PAX7','MYF5','MYOD1','MYOG','DES','CFD','APOD','ACTA2','TAGLN','MYH11','CNN1','RGS5','PDGFRB','PECAM1')
genes <- unique(c(sig$gene,markers))
stopifnot(nrow(cells)==14846L,!anyDuplicated(cells$cell_id),all(cells$all_gene_SoupX_UMI>0),
 !anyNA(match(nonzero$gene,genes)),!anyNA(match(nonzero$cell_id,cells$cell_id)),
 all(nonzero$count>0),all(nonzero$count==round(nonzero$count)))
counts <- sparseMatrix(i=match(nonzero$gene,genes),j=match(nonzero$cell_id,cells$cell_id),x=nonzero$count,
 dims=c(length(genes),nrow(cells)),dimnames=list(genes,cells$cell_id))
w <- aggregate(weight~gene,sig,mean)
libs <- c('H01045__concave','H01045__convex','H66206__concave','H66206__convex')
stopifnot(setequal(unique(cells$sample_id),libs))
result <- vector('list',nrow(sig)*length(libs)); pos <- 0L
for(i in seq_len(nrow(sig))) for(lib in libs) {
 s <- sig[i,,drop=FALSE]
 ids <- which(cells$cell_type==s$cell_type & cells$sample_id==lib)
 v <- as.numeric(counts[s$gene,ids]); n<-length(v); total<-sum(v); sorted<-sort(v,decreasing=TRUE)
 rr<-s;rr$sample_id<-lib;rr$n_cells<-n;rr$n_detected<-sum(v>0)
 rr$detection_fraction<-if(n) mean(v>0) else NA_real_
 rr$total_UMI<-total
 rr$fraction_top1<-if(total) max(v)/total else NA_real_
 rr$fraction_top5<-if(total) sum(head(sorted,5))/total else NA_real_
 rr$fraction_top1percent<-if(total) sum(head(sorted,max(1,ceiling(.01*n))))/total else NA_real_
 aa<-ann[ann$cell_type==s$cell_type,,drop=FALSE]
 rr$annotation_confidence<-paste(unique(aa$confidence),collapse=';')
 rr$reference_mapping<-paste(unique(aa$reference_mapping_evidence),collapse=';')
 rr$final_weight<-w$weight[match(s$gene,w$gene)]
 rr$detectable_GSE254300_AIS<-any(det$gene==s$gene & det$dataset=='GSE254300_AIS' & det$available)
 for(donor in c('H01045','H66206')) {
  ee<-eff$paired_logCPM_difference[eff$gene==s$gene & eff$cell_type==s$cell_type & eff$patient_id==donor]
  stopifnot(length(ee)==1L)
  rr[[paste0('donor_',donor,'_paired_logCPM_effect')]]<-ee
 }
 rr$annotation_scope<-'computational coarse label; reference mismatch does not establish misannotation'
 rr$detectable_GSE311180_AIS<-s$gene%in%other
 rr$projected_scaledTPM_full<-s$gene%in%corrected_full
 rr$projected_scaledTPM_gate<-s$gene%in%corrected_gate
 pos<-pos+1L;result[[pos]]<-rr
}
evidence<-do.call(rbind,result)
stopifnot(nrow(evidence)==600L,length(unique(evidence$gene[evidence$final_weight!=0]))==144L)
wr(evidence,'Supplementary_Table_S4_signature_evidence.tsv')
panel<-list()
for(ct in c('muscle stem/progenitor cell','smooth-muscle-like cell')) for(lib in libs) for(g in markers) {
 ids<-which(cells$cell_type==ct & cells$sample_id==lib); v<-as.numeric(counts[g,ids])
 panel[[length(panel)+1L]]<-data.frame(cell_type=ct,sample_id=lib,gene=g,n_cells=length(ids),
  n_detected=sum(v>0),fraction=mean(v>0),mean_log1p_cp10k=mean(log1p(v/cells$all_gene_SoupX_UMI[ids]*1e4)))
}
d<-do.call(rbind,panel);wr(d,'marker_panel_source.tsv')
d$gene<-factor(d$gene,levels=markers)
d$sample_id<-factor(d$sample_id,levels=libs)
p<-ggplot(d,aes(gene,sample_id,size=fraction,colour=mean_log1p_cp10k))+
 geom_point()+facet_grid(cell_type~.,scales='free_y')+
 scale_size_area(max_size=4,limits=c(0,1),name='Detected fraction')+
 scale_colour_gradient(low='#CADCE8',high='#144A70',name='Mean log1p\nCP10K')+
 guides(size=guide_legend(order=1),colour=guide_colourbar(order=2))+
 labs(x=NULL,y=NULL)+theme_bw(base_size=8,base_family='Helvetica')+
 theme(panel.grid.minor=element_blank(),axis.text.x=element_text(angle=45,hjust=1),
  strip.text.y=element_text(angle=0),legend.position='bottom')
ggsave(file.path(out,'Figure_S11_marker_evidence.pdf'),p,width=7.2,height=4.1,device=pdf)


cat('Rebuilt',nrow(evidence),'signature records and',nrow(d),'marker records from compact counts.\n')
