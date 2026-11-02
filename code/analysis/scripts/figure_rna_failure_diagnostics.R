#!/usr/bin/env Rscript
# Evidence contract: depth-dependent, directional RNA-call dropout explains the
# score deficit; post hoc sensitivities do not replace the primary pairing gates.
# Quantitative grid, R-only, 183 x 160 mm, descriptive summaries, no cell/patient
# inferential test or confidence interval; PDF/SVG/TIFF and source tables.
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(scales)})
source('scripts/lib_common.R')
rev <- 'analysis'
read <- function(name) {
  dat <- read.delim(file.path(rev,name),check.names=FALSE)
  if ('expected_metadata_pair' %in% names(dat)) {
    stopifnot(all(dat$expected_metadata_pair %in% c('True','False')))
    dat$expected_metadata_pair <- dat$expected_metadata_pair == 'True'
  }
  dat
}
loss <- read('audit_RNA_pairing_loss_decomposition.tsv')
grid <- read('audit_RNA_pairing_depth_AB_grid.tsv')
diversity <- read('audit_RNA_discordant_fragment_coordinate_diversity.tsv')
strict <- read('audit_RNA_pairing_unique_high_quality_reads.tsv')
retained <- read('audit_RNA_pairing_unique_high_quality_reads_overlap_retained.tsv')
cap <- read('audit_RNA_pairing_full_read_cap50.tsv')
stopifnot(nrow(loss)==7,nrow(grid)==1092,nrow(diversity)==14,nrow(strict)==91,nrow(retained)==91,nrow(cap) %in% c(1,6,15))
th <- theme_classic(base_size=6.5,base_family='Helvetica') + theme(
  plot.tag=element_text(face='bold',size=9),plot.title=element_text(face='bold',size=7.5),
  plot.subtitle=element_text(size=6),legend.title=element_text(size=6),legend.text=element_text(size=5.8),
  strip.background=element_blank(),strip.text=element_text(size=6.5,face='bold'),
  plot.margin=margin(4,4,4,4),legend.position='bottom',legend.key.height=unit(2,'mm'))
loss$patient_id <- factor(loss$patient_id,levels=rev(sort(loss$patient_id)))
a_data <- rbind(data.frame(patient_id=loss$patient_id,component='Het / hom',loss=loss$loss_from_het_hom),
                data.frame(patient_id=loss$patient_id,component='Opposite hom',loss=loss$loss_from_IBS0))
a <- ggplot(a_data,aes(loss,patient_id,fill=component)) + geom_col(width=.65) +
  geom_vline(xintercept=.1,linetype='dashed',colour='#AAAAAA',linewidth=.3) +
  scale_fill_manual(values=c('Het / hom'='#C56A48','Opposite hom'='#355C8A'),name=NULL) +
  scale_x_continuous(limits=c(0,.125),breaks=c(0,.05,.1),labels=label_number(accuracy=.01),expand=expansion(mult=c(0,.01))) +
  labs(tag='a',x='Contribution to (1 - relatedness)',y=NULL,title=NULL,
       subtitle=NULL) + th
b_data <- grid[grid$expected_metadata_pair & grid$patient_a %in% c('STB-084','STB-100'),]
b_data$min_het_AB <- factor(b_data$min_het_AB,levels=c(.1,.2,.3))
b_data$depth <- factor(b_data$depth,levels=c(3,7,15,30))
b <- ggplot(b_data,aes(depth,relatedness,group=min_het_AB,linetype=min_het_AB)) +
  geom_hline(yintercept=.9,linetype='dashed',colour='#AAAAAA',linewidth=.3) +
  geom_line(colour='#355C8A',linewidth=.35) + geom_point(colour='#355C8A',size=1.2) +
  geom_point(data=b_data[b_data$depth=='7' & b_data$min_het_AB=='0.2',],shape=21,size=2.8,stroke=.6,fill=NA,colour='#C56A48') +
  facet_wrap(~patient_a,nrow=1) + scale_linetype_manual(values=c('dotted','solid','longdash'),name='Minimum het AB') +
  scale_y_continuous(limits=c(.7,1.012),breaks=c(.7,.8,.9,1)) +
  labs(tag='b',x='Minimum depth in each library',y='RNA profile relatedness',title=NULL,
       subtitle=NULL) + th
c_data <- data.frame(patient_id=as.character(loss$patient_id),setting='Primary',relatedness=loss$relatedness)
for (entry in list(list(retained,'Q20 / retain'),list(strict,'Q20 / adjust'),list(cap,'Cap50'))){
  dat <- entry[[1]];dat <- dat[dat$expected_metadata_pair,]
  c_data <- rbind(c_data,data.frame(patient_id=dat$patient_a,setting=entry[[2]],relatedness=dat$relatedness))
}
c_data$setting <- factor(c_data$setting,levels=c('Primary','Q20 / retain','Q20 / adjust','Cap50'))
c_data$group <- ifelse(c_data$patient_id=='STB-084','STB-084',ifelse(c_data$patient_id=='STB-100','STB-100','Other pairs'))
c <- ggplot(c_data,aes(setting,relatedness,group=patient_id,colour=group)) +
  geom_hline(yintercept=.9,linetype='dashed',colour='#AAAAAA',linewidth=.3) +
  geom_line(linewidth=.35) + geom_point(size=1.6) +
  scale_colour_manual(values=c('Other pairs'='#BBBBBB','STB-084'='#355C8A','STB-100'='#C56A48'),name=NULL) +
  scale_y_continuous(limits=c(min(.85,min(c_data$relatedness)-.01),1.002),breaks=c(.85,.9,.95,1)) +
  labs(tag='c',x=NULL,y='RNA profile relatedness',title=NULL,
       subtitle=NULL) + th
diversity <- diversity[order(diversity$patient_id,diversity$homozygous_call_sample),]
diversity$label <- paste0(sub('^STB-','',diversity$patient_id),': ',diversity$homozygous_call_sample,'/',diversity$heterozygous_call_sample)
diversity$label <- factor(diversity$label,levels=rev(as.character(diversity$label)))
d_data <- rbind(data.frame(label=diversity$label,side='Homozygous-call side',coordinates=diversity$homo_distinct_coordinates_median,loci=diversity$original_discordant_loci),
                data.frame(label=diversity$label,side='Paired heterozygous side',coordinates=diversity$het_distinct_coordinates_median,loci=diversity$original_discordant_loci))
d <- ggplot(d_data,aes(coordinates,label,shape=side,colour=side)) +
  geom_line(aes(group=label),colour='#BBBBBB',linewidth=.25) + geom_point(size=1.6) +
  scale_shape_manual(values=c('Homozygous-call side'=16,'Paired heterozygous side'=1),name=NULL) +
  scale_colour_manual(values=c('Homozygous-call side'='#C56A48','Paired heterozygous side'='#355C8A'),name=NULL) +
  scale_x_continuous(limits=c(0,15),breaks=c(0,5,10,15),expand=expansion(mult=c(.02,.02))) +
  labs(tag='d',x='Median distinct alignment-template coordinates',y=NULL,title=NULL,
       subtitle=NULL) + th +
  theme(axis.text.y=element_text(size=6.5),legend.text=element_text(size=6))
figure <- (a+b)/(c+d)+plot_layout(widths=c(1,1.25),heights=c(1,1.1))
dir.create(file.path(rev,'figures'),recursive=TRUE,showWarnings=FALSE)
save_pub_r(figure,file.path(rev,'figures/figure_S9_RNA_failure_diagnostics'),183,160,600)
ragg::agg_png(file.path(rev,'figures/figure_S9_RNA_failure_diagnostics.png'),width=183,height=160,units='mm',res=180)
print(figure);dev.off()
write_tsv(a_data,file.path(rev,'tables/figure_S9_loss_components.tsv'))
write_tsv(b_data,file.path(rev,'tables/figure_S9_depth_AB_settings.tsv'))
write_tsv(c_data,file.path(rev,'tables/figure_S9_read_quality_cap50.tsv'))
write_tsv(d_data,file.path(rev,'tables/figure_S9_fragment_coordinate_diversity.tsv'))
