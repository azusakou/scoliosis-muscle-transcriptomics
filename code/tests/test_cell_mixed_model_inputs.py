"""Do not impute unavailable or incompletely exported cell-model inputs."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class CellMixedModelInputTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript unavailable')
    def test_real_zero_coverage_provenance_and_model_tolerance_guards(self):
        program = r'''
args <- commandArgs(trailingOnly=TRUE)
library <- function(...) stop("Input guard fixtures must not load model packages")
source(args[1]); temporary <- args[2]
md <- data.frame(cell_id=c("c1","c2","c3","c4"),cell_type="SMC",
 sample_id=c("p1c","p1v","p2c","p2v"),patient_id=c("p1","p1","p2","p2"),
 side=c("concave","convex","concave","convex"),nCount_RNA=100,nCount_SoupX=90,subsets_Mito_percent=3)
pair <- data.frame(cell_type="SMC",gene="SLN")
a <- data.frame(assay="SoupX",cell_type="SMC",gene="SLN",available=TRUE,n_cells=4L,n_nonzero=2L,
 all_gene_total_source_field="nCount_SoupX",zero_encoding="every_model_cell_explicit_including_zero")
z <- data.frame(assay="SoupX",cell_type="SMC",gene="SLN",cell_id=rev(md$cell_id),UMI=c(0,2,0,1))
reject <- function(expression) stopifnot(inherits(try(force(expression),silent=TRUE),"try-error"))
stopifnot(identical(mixed_sidecar_values(pair,md,z,a),c(1,0,2,0)))
reject(mixed_sidecar_values(pair,md,z[-1,],a))
wrong <- a; wrong$available <- FALSE; reject(mixed_sidecar_values(pair,md,z,wrong))
wrong <- z; wrong$cell_id[1] <- wrong$cell_id[2]; reject(mixed_sidecar_values(pair,md,wrong,a))
wrong <- z; wrong$gene <- "OTHER"; reject(mixed_sidecar_values(pair,md,wrong,a))
wrong <- z; wrong$UMI[1] <- .25; reject(mixed_sidecar_values(pair,md,wrong,a))
wrong <- a; wrong$all_gene_total_source_field <- "target_subset_sum"; reject(mixed_sidecar_values(pair,md,z,wrong))
mixed_validate_metadata(md)
wrong <- md; wrong$nCount_SoupX[1] <- 0; reject(mixed_validate_metadata(wrong))
wrong <- md; wrong$cell_id[1] <- wrong$cell_id[2]; reject(mixed_validate_metadata(wrong))
p <- file.path(temporary,"counts.tsv");writeLines("authentic counts",p)
binding <- list(list(path="counts.tsv",bytes=unname(file.info(p)$size),sha256=mixed_hash(p)))
mixed_check_files(temporary,binding);write("changed count",p,append=TRUE);reject(mixed_check_files(temporary,binding))
binding[[1]]$path <- "../escape";reject(mixed_check_files(temporary,binding))
mixed_close(1+1e-9,1,"existing solver tolerance")
reject(mixed_close(1+2e-7,1,"outside original criterion"))
'''
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp)/'checks.R'; script.write_text(program)
            result = subprocess.run([shutil.which('Rscript'), str(script),
                                     str(ROOT/'scripts/refit_retained_cell_mixed_models.R'), tmp],
                                    capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)


if __name__ == '__main__':
    unittest.main()
