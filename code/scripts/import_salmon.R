suppressPackageStartupMessages(library(stats))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4) stop('Usage: Rscript import_salmon.R QUANT_DIRECTORY TRANSCRIPT_GENE_MAP SAMPLE_TABLE OUTPUT.tsv.gz')
mp<-read.delim(args[2]);info<-read.delim(args[3])
res<-list()
for(id in info$sample_id){
 q<-read.delim(file.path(args[1],id,'quant.sf'));i<-match(q$Name,mp$transcript_id);k<-!is.na(i);q<-q[k,];g<-mp$gene_name[i[k]]
 x<-rowsum(cbind(TPM=q$TPM,NumReads=q$NumReads,weighted_length=q$TPM*q$EffectiveLength,length=q$EffectiveLength,n=rep(1,nrow(q))),g,reorder=TRUE)
 len<-x[,'weighted_length']/x[,'TPM'];len[!is.finite(len)]<-x[!is.finite(len),'length']/x[!is.finite(len),'n']
 res[[id]]<-data.frame(gene=rownames(x),sample_id=id,TPM=sprintf('%.17g',x[,'TPM']),NumReads=sprintf('%.17g',x[,'NumReads']),effective_length=sprintf('%.17g',len))
}
write.table(do.call(rbind,res),gzfile(args[4]),sep='\t',quote=FALSE,row.names=FALSE)
