suppressPackageStartupMessages(library(edgeR))
out <- Sys.getenv('AIS_OUTPUT', 'computed')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
rd <- function(p) read.delim(p,check.names=FALSE)
wr <- function(x,n) write.table(x,file.path(out,paste0(n,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
mat <- function(p) {d<-rd(p);x<-as.matrix(d[,-1]);rownames(x)<-d[[1]];x}
base <- ''
a <- mat(paste0(base,'results/tables/bulk/GSE254300_AIS_label_blind_projection_logCPM.tsv.gz'))
b <- mat(paste0(base,'results/tables/bulk_gse311180_requant/GSE311180_label_blind_projection_logCPM.tsv.gz'))
ia <- rd(paste0(base,'results/tables/bulk/multicellular_axis_scores.tsv'))
ib <- rd(paste0(base,'results/tables/bulk_gse311180_requant/frozen_axis_scores.tsv'))
sig <- rd('results/tables/bulk/frozen_multicellular_signature.tsv')
w <- aggregate(weight~gene,sig,mean);w<-w[w$weight!=0,];weights<-setNames(w$weight,w$gene)
exact <- function(d) {sg<-as.matrix(expand.grid(rep(list(c(-1,1)),length(d))));mean(abs(sg%*%d/length(d))>=abs(mean(d))-1e-12)}
summ <- function(d,n) {se<-sd(d)/sqrt(length(d));ci<-mean(d)+c(-1,1)*qt(.975,length(d)-1)*se;data.frame(analysis=n,n=length(d),mean=mean(d),sd=sd(d),lower=ci[1],upper=ci[2],p=exact(d))}
paired <- function(s,info) {ids<-sort(unique(info$patient_id));setNames(vapply(ids,function(id) s[info$sample_id[info$patient_id==id&info$side=='concave']]-s[info$sample_id[info$patient_id==id&info$side=='convex']],numeric(1)),ids)}
project <- function(x,genes) {z<-t(scale(t(x[genes,,drop=FALSE])));z[!is.finite(z)]<-0;colSums(z*weights[genes])/sum(abs(weights[genes]))}
shared <- intersect(intersect(rownames(a),rownames(b)),names(weights)); stopifnot(length(shared)==117)
da<-paired(project(a,shared),ia);db<-paired(project(b,shared),ib)
wr(data.frame(gene=shared,weight=weights[shared]),'shared117_weights')
wr(rbind(summ(da,'shared117_GSE254300'),summ(db,'shared117_GSE311180'),summ(c(da,db),'shared117_pooled')),'shared117_summary')
wr(data.frame(patient=c(names(da),names(db)),cohort=rep(c('GSE254300','GSE311180'),c(5,7)),difference=c(da,db)),'shared117_patient_effects')
# Reciprocal conditional tails from the original released fitted-effect grid.
grid<-rd(paste0(base,'analysis/bulk_gene_concordance_edgeR_permutation_null.tsv'))
tails<-function(g) {
 obs<-g[g[[1]]==1&g[[2]]==1,];r<-list()
 for(k in c('spearman_rho','sign_concordance')) for(j in 1:2) {
  v<-g[g[[3-j]]==1,k];ref<-obs[[k]];if(k=='sign_concordance'){v<-v-.5;ref<-ref-.5}
  r[[length(r)+1]]<-data.frame(endpoint=k,randomized_cohort=c('GSE254300','GSE311180')[j],configurations=length(v),p=mean(abs(v)>=abs(ref)-1e-12))
 };do.call(rbind,r)
}
wr(tails(grid),'original_conditional_tails')
# Effective-length-aware import: gene-summed Salmon TPMs converted to scaledTPM
# counts using each sample's total estimated fragments (tximport definition).
# No length offset is added again. This avoids dependence on isoform length.
mf<-rd('metadata/audit/gse311180_requant_manifest.tsv');mf<-mf[match(colnames(b),mf$sample_id),];stopifnot(!anyNA(mf$sample_id))
abundance<-rd('data/salmon_gene_abundance.tsv.gz')
genes<-sort(unique(abundance$gene));xs<-rawcounts<-matrix(0,length(genes),nrow(mf),dimnames=list(genes,mf$sample_id))
for(id in mf$sample_id){
 q<-abundance[abundance$sample_id==id,];stopifnot(!anyDuplicated(q$gene),all(q$TPM>=0),all(q$NumReads>=0))
 rawcounts[q$gene,id]<-q$NumReads
 xs[q$gene,id]<-q$TPM/sum(q$TPM)*sum(q$NumReads)
}
old<-mat(paste0(base,'analysis/bulk_gate_inputs/GSE311180_tested_counts.tsv.gz'))
delta<-rawcounts[rownames(old),colnames(old)]-old;wr(data.frame(sample=colnames(old),max_abs=apply(abs(delta),2,max),correlation=vapply(seq_len(ncol(old)),function(j) cor(rawcounts[rownames(old),colnames(old)[j]],old[,j]),numeric(1))),"quant_count_comparison");print(summary(as.numeric(delta)));stopifnot(min(vapply(seq_len(ncol(old)),function(j) cor(rawcounts[rownames(old),colnames(old)[j]],old[,j]),numeric(1)))>.9999999)
write.table(xs,gzfile(file.path(out,'GSE311180_scaledTPM_counts.tsv.gz')),sep='\t',quote=FALSE,col.names=NA)
ac<-mat(paste0(base,'analysis/bulk_gate_inputs/GSE254300_tested_counts.tsv.gz'))
oldcommon<-intersect(rownames(ac),rownames(old));stopifnot(length(oldcommon)==12995)
fit <- function(counts,info,universe=NULL) {
 info<-info[match(colnames(counts),info$sample_id),];design<-model.matrix(~factor(patient_id)+factor(side,levels=c('convex','concave')),info)
 y<-DGEList(counts);eligible<-rownames(counts)[filterByExpr(y,group=info$side,min.count=10)]
 if(is.null(universe)) universe<-eligible
 y<-y[universe,,keep.lib.sizes=FALSE];y<-calcNormFactors(y);y<-estimateDisp(y,design,robust=TRUE);f<-glmQLFit(y,design,robust=TRUE)
 tab<-topTags(glmQLFTest(f,coef=ncol(design)),n=Inf,sort.by='none')$table
 list(y=y,tab=tab,info=info,eligible=eligible)
}
allres<-list();uni<-list();candidates<-list();scoreout<-list()
for(mode in c('full','gate')) {
 keep<-if(mode=='full') rep(TRUE,ncol(xs)) else !mf$patient_id%in%c('STB-084','STB-100')
 xx<-xs[,keep,drop=FALSE];info<-mf[keep,];fa<-fit(ac,ia,rownames(ac));fb<-fit(xx,info,rownames(old))
 wr(data.frame(gene=rownames(fb$tab),fb$tab),paste0(mode,'_scaledTPM_edgeR'))
 native<-fb$eligible;uni[[mode]]<-data.frame(mode=mode,original=length(rownames(old)),new_eligible=length(native),lost=length(setdiff(rownames(old),native)),gained=length(setdiff(native,rownames(old))))
 wr(data.frame(gene=union(rownames(old),native),original=union(rownames(old),native)%in%rownames(old),scaledTPM=union(rownames(old),native)%in%native),paste0(mode,'_gene_eligibility'))
 rawfit<-fit(rawcounts[,keep,drop=FALSE],info,rownames(old));wr(data.frame(mode=mode,rho=cor(fa$tab[oldcommon,'logFC'],rawfit$tab[oldcommon,'logFC'],method='spearman'),agreement=mean(sign(fa$tab[oldcommon,'logFC'])==sign(rawfit$tab[oldcommon,'logFC']))),paste0(mode,'_same_quant_uncorrected_control'));wr(data.frame(gene=c('ANGPT1','LEP'),rawfit$tab[c('ANGPT1','LEP'),]),paste0(mode,'_same_quant_uncorrected_candidates'))
 # Projection has its own side-label-blind detectability filter.
 full<-DGEList(xx);full$samples$norm.factors<-fb$y$samples$norm.factors
 lc<-cpm(full,log=TRUE,prior.count=2);detect<-rowSums(cpm(full)>=1)>=max(2,ceiling(.2*ncol(xx)))
 genes<-intersect(names(weights),rownames(xx)[detect]);d<-paired(project(lc,genes),info)
 scoreout[[mode]]<-rbind(summ(d,paste0(mode,'_cohort_scaledTPM')),summ(c(paired(project(a,intersect(names(weights),rownames(a))),ia),d),paste0(mode,'_pooled_scaledTPM')))
 scoreout[[mode]]$n_signature_genes<-length(genes)
 wr(data.frame(gene=genes,weight=weights[genes]),paste0(mode,'_scaledTPM_projection_genes'))
 candidates[[mode]]<-data.frame(mode=mode,gene=c('ANGPT1','LEP'),fb$tab[c('ANGPT1','LEP'),]);candidates[[mode]]$two_candidate_BH<-p.adjust(candidates[[mode]]$PValue,'BH')
 configs<-function(n) as.matrix(expand.grid(rep(list(c(1,-1)),n)))
 permfits<-function(f) {
  ids<-unique(f$info$patient_id);sg<-configs(length(ids));ans<-matrix(NA_real_,length(oldcommon),nrow(sg))
  for(i in seq_len(nrow(sg)/2)) {
   inf<-f$info;flip<-sg[i,match(inf$patient_id,ids)]==-1;inf$side[flip]<-ifelse(inf$side[flip]=='concave','convex','concave')
   design<-model.matrix(~factor(patient_id)+factor(side,levels=c('convex','concave')),inf)
   y<-estimateDisp(f$y,design,robust=TRUE);z<-glmQLFit(y,design,robust=TRUE)
   eff<-glmQLFTest(z,coef=ncol(design))$table[oldcommon,'logFC'];ans[,i]<-eff;opp<-which(rowSums(abs(sweep(sg,2,-sg[i,],'-')))==0);ans[,opp]<- -eff
  };ans
 }
 aa<-permfits(fa);bb<-permfits(fb);gg<-expand.grid(GSE254300_configuration=seq_len(ncol(aa)),GSE311180_configuration=seq_len(ncol(bb)))
 gg$spearman_rho<-vapply(seq_len(nrow(gg)),function(i) cor(aa[,gg[i,1]],bb[,gg[i,2]],method='spearman'),numeric(1))
 gg$sign_concordance<-vapply(seq_len(nrow(gg)),function(i) mean(sign(aa[,gg[i,1]])==sign(bb[,gg[i,2]])),numeric(1))
 wr(gg,paste0(mode,'_scaledTPM_permutation_grid'));wr(tails(gg),paste0(mode,'_scaledTPM_conditional_tails'))
 allres[[mode]]<-data.frame(mode=mode,rho=gg$spearman_rho[1],agreement=gg$sign_concordance[1],rho_p=mean(abs(gg$spearman_rho)>=abs(gg$spearman_rho[1])-1e-12),agreement_p=mean(abs(gg$sign_concordance-.5)>=abs(gg$sign_concordance[1]-.5)-1e-12))
 print(allres[[mode]]);flush.console()
}
wr(do.call(rbind,uni),'scaledTPM_gene_universe');wr(do.call(rbind,candidates),'scaledTPM_candidates');wr(do.call(rbind,allres),'scaledTPM_concordance');wr(do.call(rbind,scoreout),'scaledTPM_scores')

