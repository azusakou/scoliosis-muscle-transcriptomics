"""Protect count partitions and original diagnostic decision thresholds."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class SpecificityInputTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript unavailable')
    def test_nonduplicated_partition_and_strict_decision_guards(self):
        code = r'''
args <- commandArgs(trailingOnly=TRUE)
library <- function(...) stop("Input guard fixtures must not load model packages")
source(args[1])
reject <- function(expression)stopifnot(inherits(try(force(expression),silent=TRUE),"try-error"))
genes<-paste0("g",seq_len(14967));ais<-paste0("A",seq_len(10));cs<-paste0("C",seq_len(10))
info<-data.frame(cohort_label=rep(c("AIS","CS"),each=10),processed_sample_id=c(ais,cs))
A<-matrix(1,14296,10,dimnames=list(genes[1:14296],ais));extra<-matrix(2,671,10,dimnames=list(genes[14297:14967],ais))
CS<-matrix(3,14967,10,dimnames=list(genes,cs))
joint<-specificity_reconstruct(A,extra,CS,genes,info)
stopifnot(identical(dim(joint),c(14967L,20L)),all(joint[rownames(A),ais]==1),all(joint[rownames(extra),ais]==2),all(joint[,cs]==3))
reject(specificity_reconstruct(A,extra[-1,,drop=FALSE],CS,genes,info))
wrong<-extra;rownames(wrong)[1]<-rownames(A)[1];reject(specificity_reconstruct(A,wrong,CS,genes,info))
wrong<-CS;colnames(wrong)<-rev(colnames(wrong));reject(specificity_reconstruct(A,extra,wrong,genes,info))
wrong<-extra;wrong[1,1]<-NA;reject(specificity_reconstruct(A,wrong,CS,genes,info))
x<-data.frame(logFC=1,logCPM=2,F=3,PValue=.1,FDR=.2,gene="g",contrast="fixed",n_paired_patients=5L)
specificity_compare(x,x,"same")
y<-x;y$logFC<-1+2e-7;reject(specificity_compare(y,x,"model criterion"))
y<-x;y$PValue<-.05;x$PValue<-.05-1e-12;reject(specificity_compare(y,x,"threshold decision"))
x$logFC<--1e-12;y<-x;y$logFC<-1e-12;reject(specificity_compare(y,x,"effect sign"))
y<-x;y$n_paired_patients<-4L;reject(specificity_compare(y,x,"patient count"))
'''
        with tempfile.TemporaryDirectory() as tmp:
            script=Path(tmp)/'checks.R';script.write_text(code)
            r=subprocess.run([shutil.which('Rscript'),str(script),str(ROOT/'scripts/refit_specificity_models.R')],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stdout+r.stderr)


if __name__=='__main__':unittest.main()
