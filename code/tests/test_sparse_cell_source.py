"""Reject expanded-cell corruption against actual sparse input and full totals."""
import math,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from verify_cell_audit_replay import check_sparse_row

class SparseCellSourceTests(unittest.TestCase):
    def fixture(self):
        info=dict(sample_id='s',patient_id='p',side='concave',cell_type='FAP',nCount_RNA='100',nFeature_RNA='50',subsets_Mito_percent='2.5',annotation_confidence='high')
        row=dict(info,assay='RNA',cell_id='c',gene='G',UMI='2',all_gene_total_UMI='100',cp10k='200',log1p_cp10k=str(math.log1p(200)))
        row.pop('nCount_RNA')
        return row,{'c':info},{('RNA','c','G'):2},{('RNA','G')}
    def test_source_exact_and_observed_sparse_zero(self):
        row,md,counts,available=self.fixture()
        self.assertEqual(check_sparse_row(row,md,counts,available),('RNA','c','G'))
        row.update(UMI='0',cp10k='0',log1p_cp10k='0')
        check_sparse_row(row,md,{},available)
    def test_missing_feature_cannot_become_zero(self):
        row,md,counts,available=self.fixture();row.update(UMI='0',cp10k='0',log1p_cp10k='0')
        with self.assertRaisesRegex(ValueError,'Unavailable'):check_sparse_row(row,md,{},set())
    def test_count_label_full_total_and_normalization_corruption_rejected(self):
        row,md,counts,available=self.fixture()
        for field,value in [('UMI','3'),('UMI','2.00000000001'),('UMI','nan'),('patient_id','wrong'),('cell_id','unknown'),('all_gene_total_UMI','99'),('cp10k','201'),('log1p_cp10k','0'),('nFeature_RNA','50.0000000001'),('subsets_Mito_percent','3'),('annotation_confidence','low')]:
            with self.subTest(field=field,value=value),self.assertRaises(ValueError):
                check_sparse_row(dict(row,**{field:value}),md,counts,available)
if __name__=='__main__':unittest.main()
