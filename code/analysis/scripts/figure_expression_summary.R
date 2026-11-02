#!/usr/bin/env Rscript
# Contract: show limited patient-level evidence, discovery heterogeneity and
# exact null discreteness without implying multicellular coordination or equivalence.
# R exclusively; 183-mm editable PDF/SVG plus 600-dpi TIFF and R-rendered PNG.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(ggrepel); library(scales)})
source("scripts/lib_common.R")
out <- "analysis/figures"
tab <- "analysis/tables"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(tab, recursive=TRUE, showWarnings=FALSE)
rd <- function(p) {p <- analysis_input(p); read.delim(if(grepl("\\.gz$",p)) gzfile(p) else p, check.names=FALSE)}
save_figure <- function(p, name, h=130) {
  save_pub_r(p, file.path(out,name),183,h,600)
  ragg::agg_png(file.path(out,paste0(name,".png")),width=183,height=h,units="mm",res=180)
  print(p); dev.off()
}
col <- c(convex="#355C8A",concave="#C56A48",D1="#355C8A",D2="#C56A48",CS="#C56A48",AIS="#355C8A")
short <- c("fibro-adipogenic progenitor"="FAP","muscle stem/progenitor cell"="MuSC",
 "macrophage/monocyte"="Mac/Mono","T/NK cell"="T/NK","B/plasma cell"="B/plasma",
 "endothelial cell"="Endothelial","pericyte"="Pericyte","smooth-muscle-like cell"="SMC-like",
 "Schwann/perineurial cell"="Schwann/perineurial")
cellcol <- setNames(c("#D7A554","#75649A","#B85450","#587FA7","#8FB2CC",
 "#3E9891","#C88439","#8B7351","#8A8E93"),names(short))
th <- function(s=7) theme_classic(base_size=s,base_family="Helvetica") + theme(
 plot.title=element_text(face="bold",size=s+0.5),plot.subtitle=element_text(size=s-0.5),
 plot.tag=element_text(face="bold",size=9),legend.title=element_text(size=s-0.3),
 legend.text=element_text(size=s-0.7),strip.background=element_blank(),
 strip.text=element_text(face="bold"),plot.margin=margin(5,6,5,6))

# Atlas: original embedding, explicit donor/side diagnostics and lineage legend.
u <- rd("results/tables/annotation/annotated_umap.tsv")
u$donor <- factor(u$patient_id,levels=c("H01045","H66206"),labels=c("D1","D2"))
u$side <- factor(u$side,levels=c("convex","concave"),labels=c("Convex","Concave"))
ctr <- aggregate(cbind(umapunintegrated_1,umapunintegrated_2)~cell_type,u,median)
ctr$label <- unname(short[ctr$cell_type])
fapctr <- aggregate(cbind(umapunintegrated_1,umapunintegrated_2)~cell_type+donor,u[u$cell_type=="fibro-adipogenic progenitor",],median)
fapctr$label <- paste0("FAP (",fapctr$donor,")")
ctr <- rbind(ctr[ctr$cell_type!="fibro-adipogenic progenitor",],fapctr[,names(ctr)])
a <- ggplot(u,aes(umapunintegrated_1,umapunintegrated_2,colour=cell_type)) +
 geom_point(size=.1,alpha=.55) + geom_label_repel(data=ctr,aes(label=label),size=1.9,
 colour="#252A31",fill="white",label.size=0,seed=722100,max.overlaps=Inf,box.padding=.2,min.segment.length=0) +
 scale_colour_manual(values=cellcol,labels=short,guide="none") + coord_equal() +
 labs(tag="a",x="UMAP 1",y="UMAP 2",title=NULL,subtitle=NULL) + th()
b <- ggplot(u,aes(umapunintegrated_1,umapunintegrated_2,colour=donor)) +
 geom_point(size=.1,alpha=.55) + facet_wrap(~side,nrow=1) +
 scale_colour_manual(values=col[c("D1","D2")],name="Donor",labels=c("D1: H01045","D2: H66206")) +
 coord_equal() + labs(tag="b",x="UMAP 1",y="UMAP 2",title=NULL,
 subtitle=NULL) + th(6.7) + theme(legend.position="bottom") +
 guides(colour=guide_legend(override.aes=list(size=2,alpha=1)))
cnt <- rd("results/tables/composition/cell_type_counts_and_proportions.tsv")
cnt$cell_type <- factor(cnt$cell_type,levels=names(short))
cnt$sample <- factor(paste(cnt$patient_id,cnt$side,sep="\n"),levels=c("H01045\nconvex","H01045\nconcave","H66206\nconvex","H66206\nconcave"))
c <- ggplot(cnt,aes(sample,proportion,fill=cell_type)) + geom_col(width=.75) +
 scale_fill_manual(values=cellcol,labels=short,name="Lineage",drop=FALSE) +
 scale_y_continuous(labels=percent_format(accuracy=1),expand=c(0,0)) +
 labs(tag="c",x=NULL,y="Proportion of retained cells",title=NULL,
 subtitle=NULL) + th(6.6) +
 theme(legend.position="bottom",axis.text.x=element_text(size=5.8)) + guides(fill=guide_legend(ncol=3,byrow=TRUE))
ce <- rd("results/tables/composition/paired_composition_effects.tsv")
ce$label <- factor(unname(short[ce$cell_type]),levels=rev(unname(short)))
ce$donor <- factor(ce$patient_id,levels=c("H01045","H66206"),labels=c("D1","D2"))
ce$effect_pp <- 100*ce$concave_minus_convex
ce$missing_label <- ifelse(is.finite(ce$effect_pp),"","not estimable")
d <- ggplot(ce[is.finite(ce$effect_pp),],aes(effect_pp,label,colour=donor,shape=donor)) +
 geom_vline(xintercept=0,linewidth=.3) + geom_point(position=position_dodge(width=.45),size=1.9) +
 geom_text(data=ce[!is.finite(ce$effect_pp),],aes(x=-4,label=missing_label),size=1.7,hjust=0,
 position=position_dodge(width=.45),show.legend=FALSE) +
 scale_colour_manual(values=col[c("D1","D2")],name="Donor") + scale_shape_manual(values=c(16,17),name="Donor") +
 labs(tag="d",x="Concave minus convex (percentage points)",y=NULL,title=NULL,
 subtitle=NULL) + th(6.6) + theme(legend.position="bottom")
save_figure((a|b)/(c|d),"figure_atlas",174)
write_tsv(cnt,file.path(tab,"figure_atlas_absolute_composition.tsv"))
write_tsv(ce,file.path(tab,"figure_atlas_paired_composition.tsv"))

# Programmes: numerical effects and a heuristic 50% independent-sign reference.
prog <- rd("results/figures/source_data/Figure_3_donor_pathways.tsv")
pathway_labels <- c(HALLMARK_ANGIOGENESIS="Angiogenesis",
 HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION="EMT",
 HALLMARK_ESTROGEN_RESPONSE_EARLY="Estrogen response (early)",
 HALLMARK_HYPOXIA="Hypoxia",HALLMARK_IL6_JAK_STAT3_SIGNALING="IL6-JAK-STAT3",
 HALLMARK_INFLAMMATORY_RESPONSE="Inflammatory response",HALLMARK_MYOGENESIS="Myogenesis",
 HALLMARK_NOTCH_SIGNALING="Notch signalling",HALLMARK_OXIDATIVE_PHOSPHORYLATION="Oxidative phosphorylation",
 HALLMARK_TGF_BETA_SIGNALING="TGF-beta signalling",HALLMARK_TNFA_SIGNALING_VIA_NFKB="TNF-alpha via NF-kappaB",
 HALLMARK_WNT_BETA_CATENIN_SIGNALING="WNT/beta-catenin")
stopifnot(setequal(unique(prog$gene_set),names(pathway_labels)),nrow(prog)==120L)
prog$pathway <- unname(pathway_labels[prog$gene_set])
prog$donor <- factor(prog$patient_id,levels=c("H01045","H66206"),labels=c("D1","D2"))
prog$cell_short <- factor(prog$cell_short,levels=c("Endothelial","FAP","Mac/Mono","MuSC","SMC-like"))
path_levels <- unique(prog$pathway)
prog$pathway <- factor(prog$pathway,levels=rev(path_levels))
mx <- max(abs(prog$paired_difference))
a <- ggplot(prog,aes(donor,pathway,fill=paired_difference)) + geom_tile(colour="white",linewidth=.4) +
 geom_text(aes(label=sprintf("%.2f",paired_difference)),size=1.75) + facet_grid(.~cell_short) +
 scale_fill_gradient2(low=col["convex"],high=col["concave"],mid="#F7F7F5",limits=c(-mx,mx),name="Mean gene-wise z-score\nconcave minus convex") +
 labs(tag="a",x="Discovery donor",y=NULL,title=NULL,
 subtitle=NULL) + th(6.6) +
 theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom")
conc <- rd("results/figures/source_data/Figure_3_concordance.tsv")
conc$n_concordant <- round(conc$concordant*12)
b <- ggplot(conc,aes(concordant_pct,reorder(cell_short,concordant_pct))) +
 geom_vline(xintercept=50,linetype=2,colour="#666666",linewidth=.4) +
 geom_point(size=2.4,colour="#3A8F8A") + geom_text(aes(label=paste0(n_concordant,"/12")),hjust=-.4,size=2) +
 scale_x_continuous(limits=c(0,105),breaks=c(0,50,100),labels=function(x)paste0(x,"%")) +
 labs(tag="b",x="Same sign in both donors",y=NULL,title=NULL,
 subtitle=NULL) + th(6.5)
save_figure(a|b+plot_layout(widths=c(1.8,1)),"figure_programmes",126)
write_tsv(prog,file.path(tab,"figure_programmes_source.tsv"));write_tsv(conc,file.path(tab,"figure_programmes_concordance.tsv"))

# Signature construction: selected extreme tail, with actual average abundance.
sig <- rd("results/tables/bulk/frozen_multicellular_signature.tsv")
conf <- rd("results/tables/bulk/frozen_multicellular_signature_direction_conflicts.tsv")
sig$cell_short <- unname(short[sig$cell_type]); sig$direction <- ifelse(sig$weight>0,"Concave-high","Convex-high")
sn <- as.data.frame(table(sig$cell_short,sig$direction));names(sn)<-c("lineage","direction","n")
sn$plot_n<-ifelse(sn$direction=="Concave-high",sn$n,-sn$n)
a <- ggplot(sn,aes(plot_n,lineage,fill=direction))+geom_vline(xintercept=0,linewidth=.3)+geom_col(width=.65)+
 scale_fill_manual(values=c("Concave-high"=unname(col["concave"]),"Convex-high"=unname(col["convex"])),name=NULL)+scale_x_continuous(labels=abs)+
 labs(tag="a",x="Selected lineage-gene entries",y=NULL,title=NULL,subtitle=NULL)+th(6.7)+theme(legend.position="bottom")
st <- sig[!sig$gene%in%conf$gene,]
st <- do.call(rbind,lapply(split(st,st$cell_type),function(z)head(z[order(abs(z$mean_donor_effect),decreasing=TRUE),],2)))
de <- rd("results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz")
st <- merge(st,de[,c("gene","cell_type","logCPM")],by=c("gene","cell_type"))
st$label <- sprintf("%s · %s (%.1f)",st$gene,st$cell_short,st$logCPM)
b <- ggplot(st,aes(mean_donor_effect,reorder(label,mean_donor_effect)))+geom_vline(xintercept=0,linewidth=.3)+
 geom_segment(aes(x=0,xend=mean_donor_effect,yend=reorder(label,mean_donor_effect)),colour="#DDE1E5",linewidth=.5)+
 geom_point(colour="#355C8A",size=1.9)+labs(tag="b",x="Mean donor log2 CPM difference",y=NULL,
 title=NULL,subtitle=NULL)+th(6.3)
sc <- rd("results/tables/bulk/multicellular_axis_scores.tsv");sc$side<-factor(sc$side,levels=c("convex","concave"))
mn <- aggregate(score~side,sc,mean)
c <- ggplot(sc,aes(side,score,group=patient_id))+geom_hline(yintercept=0,colour="#DDE1E5",linewidth=.3)+
 geom_line(colour="#90969C",linewidth=.45)+geom_point(aes(shape=side),size=2)+
 geom_point(data=mn,aes(side,score,shape="Mean"),inherit.aes=FALSE,size=2.8)+
 scale_shape_manual(values=c(convex=1,concave=16,Mean=18),name=NULL,labels=c(convex="Convex",concave="Concave",Mean="Mean"))+
 labs(tag="c",x=NULL,y="Mean signed gene-wise z-score",title=NULL,subtitle=NULL)+th(6.7)+theme(legend.position="bottom")
null <- rd("analysis/bulk_GSE254300_corrected_random_signature_null.tsv")
obs <- mean(rd("results/tables/bulk/multicellular_axis_paired_effects.tsv")$paired_difference)
d <- ggplot(null,aes(random_mean))+geom_histogram(bins=35,fill="#DDE1E5",colour="white")+
 geom_vline(xintercept=obs,colour="#B23A48",linewidth=.7)+labs(tag="d",x="Mean paired difference",y="Random signatures",
 title=NULL,subtitle=NULL)+th(6.5)
save_figure((a|c)/(b|d),"figure_score",158)
write_tsv(st,file.path(tab,"figure_score_selected_genes_abundance.tsv"))

# Validation: fixed shared AIS/CS scales, named symbols, patient-level uncertainty.
rep <- rd("results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz")
edgeR_null <- rd("analysis/bulk_gene_concordance_edgeR_permutation_summary.tsv")
a <- ggplot(rep,aes(logFC.GSE254300,logFC.GSE311180_requant))+geom_hline(yintercept=0,colour="#DDE1E5",linewidth=.3)+
 geom_vline(xintercept=0,colour="#DDE1E5",linewidth=.3)+geom_bin_2d(bins=45)+
 scale_fill_gradient(low="#EEF1F4",high="#294B73",trans="sqrt",guide="none")+coord_fixed()+
 labs(tag="a",x="GSE254300 fitted edgeR log2 fold change",y="GSE311180 fitted edgeR log2 fold change",title=NULL,
 subtitle=NULL)+th(6.5)+
 theme(plot.subtitle=element_text(size=6,lineheight=1.05),plot.margin=margin(5,12,5,6))
p12<-rd("results/tables/meta_validation_requant/frozen_axis_12_patient_effects.tsv")
p12$cohort<-ifelse(p12$cohort=="GSE254300","GSE254300","GSE311180")
rna_pairs<-rd("analysis/audit_gse311180_rna_expected_pairs.tsv")
failed<-rna_pairs$patient_id[rna_pairs$interpretation=="does_not_meet_consistency_gates"]
stopifnot(setequal(failed,c("STB-084","STB-100")))
p12$rna_gate_status<-ifelse(p12$cohort=="GSE254300","not_assessed_here",ifelse(p12$patient_id%in%failed,"nonpassing","passing"))
p12$display_label<-paste0(p12$patient_id,ifelse(p12$rna_gate_status=="nonpassing","*",""))
p12$patient<-factor(p12$display_label,levels=p12$display_label[order(p12$paired_difference)])
b<-ggplot(p12,aes(patient,paired_difference,colour=cohort))+geom_hline(yintercept=0,linewidth=.3)+
 geom_segment(aes(xend=patient,y=0,yend=paired_difference),linewidth=.6)+geom_point(size=1.9)+
 scale_colour_manual(values=c(GSE254300=unname(col["convex"]),GSE311180="#7E9DBF"),name="Cohort")+
 labs(tag="b",x=NULL,y="Mean signed z-score: concave minus convex",title=NULL,
 subtitle=NULL,
 caption="* Below recorded RNA relatedness gate; original scores retained.")+th(6.5)+
 theme(axis.text.x=element_text(angle=50,hjust=1),legend.position="bottom")
cs<-rd("results/tables/disease_control_gse254300_cs/frozen_axis_joint_scores.tsv")
cs$side<-factor(cs$side,levels=c("convex","concave"));cs$cohort<-factor(cs$cohort,levels=c("AIS","CS"))
cm<-aggregate(score~cohort+side,cs,mean)
c<-ggplot(cs,aes(side,score,group=pair_id,colour=cohort))+geom_hline(yintercept=0,colour="#DDE1E5",linewidth=.3)+
 geom_line(linewidth=.4,alpha=.5)+geom_point(aes(shape=side),size=1.8)+
 geom_point(data=cm,aes(side,score,colour=cohort,shape="Mean"),inherit.aes=FALSE,size=2.7)+
 facet_grid(.~cohort,scales="fixed")+scale_colour_manual(values=col[c("AIS","CS")],guide="none")+
 scale_shape_manual(values=c(convex=1,concave=16,Mean=18),name=NULL,labels=c(convex="Convex",concave="Concave",Mean="Mean"))+
 labs(tag="c",x=NULL,y="Mean signed gene-wise z-score",title=NULL,
 subtitle=NULL,
 caption="Interaction P = 0.762; failure to detect an interaction does not establish equivalence.")+th(6.5)+theme(legend.position="bottom")
mo<-rd("results/tables/disease_control_gse254300_cs/frozen_axis_model_sensitivity.tsv")
mo$lower<-mo$estimate_concave_minus_convex-qt(.975,mo$residual_df)*mo$se
mo$upper<-mo$estimate_concave_minus_convex+qt(.975,mo$residual_df)*mo$se
mo$model<-factor(mo$model,levels=c("biology_plus_composition_PC1","biology_plus_library_QC","biology_only"),labels=c("+ composition PC1","+ library size","Patient + side"))
d<-ggplot(mo,aes(estimate_concave_minus_convex,model,colour=cohort))+geom_vline(xintercept=0,linewidth=.3)+
 geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=0,position=position_dodge(width=.35),linewidth=.5)+
 geom_point(position=position_dodge(width=.35),size=1.8)+scale_colour_manual(values=col[c("AIS","CS")],name="Diagnosis")+
 labs(tag="d",x="Adjusted paired effect (95% t confidence interval)",y=NULL,title=NULL,
 subtitle=NULL)+th(6.3)+theme(legend.position="bottom")
save_figure((a|b)/(c|d)+plot_layout(widths=c(1,1)),"figure_validation",166)
write_tsv(mo,file.path(tab,"figure_validation_t_confidence_intervals.tsv"))
write_tsv(rep,file.path(tab,"figure_validation_edgeR_gene_effects.tsv"))
write_tsv(p12[,c("cohort","patient_id","paired_difference","rna_gate_status","display_label")],
 file.path(tab,"figure_validation_original_patient_scope.tsv"))

# Robustness: exact configurations are discrete; separate R2 ranges.
perm<-rd("results/tables/robustness/bulk_axis_permutation_distribution.tsv")
perm<-perm[order(perm$permuted_mean),];perm$rank<-seq_len(nrow(perm))
a<-ggplot(perm,aes(permuted_mean,0))+geom_point(size=1.6,colour="#7A818A",position=position_jitter(height=.08,width=0,seed=254300))+
 geom_vline(xintercept=obs,colour="#B23A48",linewidth=.6)+scale_y_continuous(breaks=NULL)+
 labs(tag="a",x="Mean paired effect under exact side flips",y=NULL,title=NULL,
 subtitle=NULL)+th(6.5)
poolednull<-rd("analysis/bulk_pooled_matched_random_signature_null.tsv")
pooledobs<-mean(p12$paired_difference)
b<-ggplot(poolednull,aes(pooled12_mean))+geom_histogram(bins=40,fill="#DDE1E5",colour="white")+
 geom_vline(xintercept=pooledobs,colour="#B23A48",linewidth=.7)+labs(tag="b",x="Pooled matched random-signature mean effect",y="Random signatures",
 title=NULL,subtitle=NULL)+th(6.5)
loo<-rd("analysis/bulk_pooled_leave_one_patient_out.tsv")
loo$omitted_patient<-factor(loo$omitted_patient,levels=rev(loo$omitted_patient))
c<-ggplot(loo,aes(mean_paired_difference,omitted_patient))+geom_vline(xintercept=0,linewidth=.3)+
 geom_vline(xintercept=pooledobs,linetype=2,colour="#B23A48",linewidth=.4)+
 geom_errorbar(aes(xmin=CI95_lower,xmax=CI95_upper),orientation="y",width=0,linewidth=.4,colour="#7A818A")+
 geom_point(size=1.6,colour="#3A8F8A")+
 labs(tag="c",x="Mean retained effect (95% t confidence interval)",y="Omitted patient",title=NULL,subtitle=NULL)+th(6.3)
mod<-rd("results/tables/bulk/multicellular_axis_model_sensitivity.tsv")
mod<-mod[is.finite(mod$estimate_concave_minus_convex),]
if(!"residual_df"%in%names(mod))mod$residual_df<-c(4,3,3)[match(mod$model,c("biology_only","biology_plus_library_QC","biology_plus_composition_PC1"))]
mod$lower<-mod$estimate_concave_minus_convex-qt(.975,mod$residual_df)*mod$se
mod$upper<-mod$estimate_concave_minus_convex+qt(.975,mod$residual_df)*mod$se
mod$label<-factor(mod$model,levels=c("biology_plus_composition_PC1","biology_plus_library_QC","biology_only"),labels=c("+ composition PC1","+ library size","Patient + side"))
d<-ggplot(mod,aes(estimate_concave_minus_convex,label))+geom_vline(xintercept=0,linewidth=.3)+
 geom_errorbar(aes(xmin=lower,xmax=upper),orientation="y",width=0,linewidth=.5)+geom_point(size=2,colour=col["convex"])+
 labs(tag="d",x="Paired effect (95% t confidence interval)",y=NULL,title=NULL)+th(6.3)
integ<-rd("results/tables/robustness/integration_embedding_sensitivity_summary.tsv")
il<-rbind(data.frame(method=integ$method,metric="Sample R2",value=integ$sample_r2),data.frame(method=integ$method,metric="Side variance fraction",value=integ$side_semipartial_r2))
e<-ggplot(il,aes(value,method))+geom_point(size=2,colour="#355C8A")+facet_wrap(~metric,scales="free_x")+
 labs(tag="e",x=expression(R^2),y=NULL,title=NULL,
 subtitle=NULL)+th(6.4)
save_figure((a|b)/(c|d)/e,"figure_robustness",164)
write_tsv(perm,file.path(tab,"figure_exact_permutation_configurations.tsv"))
write_tsv(mod,file.path(tab,"figure_robustness_t_confidence_intervals.tsv"))

# Supplement S5: separate method scales, remove misleading combined-effect units.
tf<-rd("results/tables/programs/scrna_gtrd_tf_target_scores_top_variable.tsv")
te<-merge(tf[tf$side=="concave",],tf[tf$side=="convex",],by=c("cell_type","tf_target_set","patient_id"),suffixes=c(".cc",".cv"))
te$effect<-te$score.cc-te$score.cv
selected<-rd("results/figures/source_data/Figure_5_TF.tsv")
te<-merge(te,selected[,c("cell_type","tf_target_set")],by=c("cell_type","tf_target_set"))
te$label<-paste(sub("_TARGET_GENES$","",te$tf_target_set),unname(short[te$cell_type]),sep=" · ")
te$donor<-factor(te$patient_id,levels=c("H01045","H66206"),labels=c("D1","D2"))
a<-ggplot(te,aes(effect,label,colour=donor,shape=donor))+geom_vline(xintercept=0,linewidth=.3)+geom_point(size=2,position=position_dodge(width=.4))+
 scale_colour_manual(values=col[c("D1","D2")],name="Donor")+scale_shape_manual(values=c(16,17),name="Donor")+
 labs(tag="a",x="Mean target-gene z-score difference",y=NULL,title=NULL,
 subtitle=NULL)+th(6.5)+theme(legend.position="bottom")
lr<-rd("results/tables/communication/cross_method_donor_aware_edges.tsv"); keep<-lr[lr$retained,]
ln<-aggregate(retained~source+target,keep,length);names(ln)[3]<-"n"
ln$source<-unname(short[ln$source]);ln$target<-unname(short[ln$target])
b<-ggplot(ln,aes(source,target,size=n))+geom_point(colour="#355C8A",alpha=.8)+scale_size_area(max_size=9,name="Retained edges")+
 labs(tag="b",x="Sender",y="Receiver",title=NULL,subtitle=NULL)+
 th(6.5)+theme(axis.text.x=element_text(angle=90,hjust=1),legend.position="bottom")
top<-rd("results/figures/source_data/Figure_5_communication.tsv")
top$label<-paste0(top$ligand,"-",top$receptor,"  ",unname(short[top$source])," > ",unname(short[top$target]))
tl<-rbind(data.frame(label=top$label,method="NATMI",effect=top$mean.natmi),data.frame(label=top$label,method="CellPhoneDB",effect=top$mean.cellphonedb))
c<-ggplot(tl,aes(effect,label))+geom_vline(xintercept=0,linewidth=.3)+geom_point(size=1.9,colour="#355C8A")+
 facet_wrap(~method,scales="free_x",labeller=as_labeller(c(CellPhoneDB="CellPhoneDB\ninteraction-expression score",NATMI="NATMI\nedge-specificity score")))+
 scale_x_continuous(breaks=function(x)pretty(x,n=4),labels=label_number(accuracy=.01))+
 labs(tag="c",x="Concave minus convex (native score units)",y=NULL,title=NULL,
 subtitle=NULL)+th(6.1)+
 theme(panel.spacing.x=unit(12,"mm"),strip.text=element_text(lineheight=1.05))
save_figure((a|b)/c,"figure_S5_communication",171)
write_tsv(te,file.path(tab,"figure_S5_target_set_donor_effects.tsv"));write_tsv(tl,file.path(tab,"figure_S5_separate_method_effects.tsv"))
write_tsv(data.frame(candidate_rows=nrow(lr),retained_rows=nrow(keep)),file.path(tab,"figure_S5_candidate_space.tsv"))

# Supplement S6: deposited and requantified units remain separate.
dep<-rd("metadata/audit/gse311180_expression_identity_diagnostics.tsv")
raw<-rd("results/tables/bulk_gse311180_requant/GSE311180_tissue_identity.tsv")
dep<-dep[dep$sample_id%in%raw$sample_id,]
dep$sample<-factor(dep$sample_id,levels=raw$sample_id);raw$sample<-factor(raw$sample_id,levels=raw$sample_id)
dep_l<-rbind(data.frame(sample=dep$sample,side=dep$side,marker="CKM",value=log2(dep$CKM+1)),data.frame(sample=dep$sample,side=dep$side,marker="ACAN",value=log2(dep$ACAN+1)))
raw_l<-rbind(data.frame(sample=raw$sample,side=raw$side,marker="CKM",value=raw$CKM),data.frame(sample=raw$sample,side=raw$side,marker="ACAN",value=raw$ACAN))
a<-ggplot(dep_l,aes(sample,value,colour=marker,shape=side))+geom_point(size=2)+scale_colour_manual(values=c(CKM="#355C8A",ACAN="#C56A48"),name="Tissue marker")+
 scale_shape_manual(values=c(concave=16,convex=1),name="Deposited side")+
 labs(tag="a",x=NULL,y="log2(deposited TPM + 1)",title=NULL,
 subtitle=NULL)+th(6.5)+theme(axis.text.x=element_text(angle=55,hjust=1),legend.position="bottom")
b<-ggplot(raw_l,aes(sample,value,colour=marker,shape=side))+geom_point(size=2)+scale_colour_manual(values=c(CKM="#355C8A",ACAN="#C56A48"),name="Tissue marker")+
 scale_shape_manual(values=c(concave=16,convex=1),name="Metadata side")+
 labs(tag="b",x="Library (same order in both panels)",y="TMM-normalised log2 CPM",title=NULL,
 subtitle=NULL)+th(6.5)+theme(axis.text.x=element_text(angle=55,hjust=1),legend.position="bottom")
save_figure(a/b,"figure_S6_deposit_identity",137)
write_tsv(dep_l,file.path(tab,"figure_S6_deposited_identity.tsv"));write_tsv(raw_l,file.path(tab,"figure_S6_requantified_identity.tsv"))
message("Expression summary figures rendered from retained source data.")
