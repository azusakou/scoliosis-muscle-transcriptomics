import csv,hashlib,importlib,sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
r=importlib.import_module('retained_reproduction')

class RetainedReproductionTests(unittest.TestCase):
    def table(self,d,name,text):
        p=Path(d)/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text);return p
    def test_observed_last_digit_difference_preserves_rows_text_and_counts(self):
        with tempfile.TemporaryDirectory() as d:
            a=self.table(d,'a.tsv','id\tn\tvalue\nP\t16469\t0.8941233608547838\n')
            b=self.table(d,'b.tsv','id\tn\tvalue\nP\t16469\t0.894123360854784\n')
            self.assertEqual(r.compare_tsv(a,b)['numeric_cells_differing'],1)
    def test_changed_values_text_structure_and_large_integer_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            original='id\tn\tvalue\nP\t100000000000000000000\t0.9\n'
            a=self.table(d,'a.tsv',original)
            for text in [original.replace('0.9','0.9001'),original.replace('P\t','X\t'),original+'P\t1\t0.9\n',original.replace('value','other'),original.replace('100000000000000000000','100000000000000000001')]:
                b=self.table(d,'b.tsv',text)
                with self.assertRaises(ValueError):r.compare_tsv(a,b)
    def test_staging_copies_only_manifest_and_rejects_changed_reference(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=self.table(d,'code.py','print(1)\n');self.table(d,'private.txt','not distributed')
            with (root/'SHA256SUMS.tsv').open('w',newline='') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256']);w.writerow(['code.py',p.stat().st_size,r.digest(p)])
            staged=r.stage(root);self.assertEqual((staged/'code.py').read_bytes(),p.read_bytes());self.assertFalse((staged/'private.txt').exists())
            p.write_text('print(2)\n')
            with self.assertRaises(ValueError):r.stage(root)
    def test_comparison_never_tolerates_changed_immutable_code(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=self.table(d,'code.py','print(1)\n')
            with (root/'SHA256SUMS.tsv').open('w',newline='') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256']);w.writerow(['code.py',p.stat().st_size,r.digest(p)])
            staged=r.stage(root);(staged/'code.py').write_text('print(2)\n')
            with self.assertRaises(ValueError):r.compare(root,staged)
            self.assertEqual(p.read_text(),'print(1)\n')
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=self.table(d,'analysis/audit_RNA.tsv','r\n0.9\n')
            with (root/'SHA256SUMS.tsv').open('w',newline='') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256']);w.writerow(['analysis/audit_RNA.tsv',p.stat().st_size,r.digest(p)])
            staged=r.stage(root);(staged/'analysis/audit_RNA.tsv').write_text('r\n0.9000000000000001\n')
            with self.assertRaises(ValueError):r.compare(root,staged)

    def test_named_derived_expression_audits_allow_only_last_bits(self):
        for name in sorted(r.DERIVED_EXPRESSION_AUDITS):
            with tempfile.TemporaryDirectory() as d:
                root=Path(d);original='sample\tshared_genes\trho\nDR73\t22816\t0.4560615816462129\n'
                p=self.table(d,name,original)
                with (root/'SHA256SUMS.tsv').open('w',newline='') as f:
                    w=csv.writer(f,delimiter='\t');w.writerow(['path','bytes','sha256']);w.writerow([name,p.stat().st_size,r.digest(p)])
                staged=r.stage(root);target=staged/name
                target.write_text(original.replace('0.4560615816462129','0.456061581646213'))
                self.assertEqual(len(r.compare(root,staged)['numerically_equivalent_tables']),1)
                for bad in [original.replace('22816','22817'),original.replace('DR73','DR74'),original.replace('0.4560615816462129','0.4561')]:
                    target.write_text(bad)
                    with self.assertRaises(ValueError):r.compare(root,staged)
                self.assertEqual(p.read_text(),original)
    def test_explicit_model_fields_accept_existing_solver_tolerance_only(self):
        with tempfile.TemporaryDirectory() as d:
            original='gene\tn_pairs\tlogFC\tFDR\tF\nANGPT1\t5\t0.72\t0.0623138289\t0\n'
            a=self.table(d,'a.tsv',original)
            b=self.table(d,'b.tsv',original.replace('0.0623138289','0.0623138339').replace('\t0\n','\t0.000000005\n'))
            columns={'logFC','FDR','F'}
            self.assertEqual(r.compare_tsv(a,b,columns)['numeric_cells_differing'],2)
            # The default arithmetic comparison is deliberately not relaxed.
            with self.assertRaises(ValueError):r.compare_tsv(a,b)
            for text in [original.replace('0.72','0.72001'),original.replace('\t5\t','\t6\t'),original.replace('ANGPT1','LEP')]:
                b=self.table(d,'b.tsv',text)
                with self.assertRaises(ValueError):r.compare_tsv(a,b,columns)
    def test_model_tolerance_never_accepts_changed_significance_or_effect_sign(self):
        with tempfile.TemporaryDirectory() as d:
            for column,old,new in [('FDR','0.04999999','0.05000001'),('logFC','0.00000001','-0.00000001'),('logFC','0','-0.00000001'),('QL95_lower','-0.00000001','0.00000001'),('PValue','1','1.00000001')]:
                a=self.table(d,'a.tsv',f'gene\t{column}\nG\t{old}\n')
                b=self.table(d,'b.tsv',f'gene\t{column}\nG\t{new}\n')
                with self.assertRaises(ValueError):r.compare_tsv(a,b,{column})
    def test_model_allowlist_is_confined_to_five_named_outputs(self):
        self.assertEqual(len(r.GATE_MODEL_COLUMNS),5)
        for columns in r.GATE_MODEL_COLUMNS.values():
            self.assertTrue({'gene','n_pairs','positive_pairs','normalization',
                'patient_paired_logCPM_mean','patient_paired_logCPM_SD',
                'patient_t_CI95_lower','patient_t_CI95_upper'}.isdisjoint(columns))
    def test_arithmetic_tolerance_never_changes_significance_decision(self):
        with tempfile.TemporaryDirectory() as d:
            for column in ['p_value','exact_p','FDR','BH_two_fixed_candidates_within_normalization',
                           'PValue.GSE254300','FDR.AIS','p_val_adj','adj.P.Val','Q_p','qvalue']:
                a=self.table(d,'a.tsv',f'gene\t{column}\nG\t0.0499999999999\n')
                b=self.table(d,'b.tsv',f'gene\t{column}\nG\t0.0500000000001\n')
                with self.assertRaises(ValueError):r.compare_tsv(a,b)
    def test_arithmetic_tolerance_never_changes_patient_interval_sign(self):
        with tempfile.TemporaryDirectory() as d:
            for column,old,new in [('patient_t_CI95_lower','-1e-14','1e-14'),
                ('patient_t_CI95_upper','0','1e-14'),('patient_paired_logCPM_mean','0','-1e-14')]:
                a=self.table(d,'a.tsv',f'gene\t{column}\nG\t{old}\n')
                b=self.table(d,'b.tsv',f'gene\t{column}\nG\t{new}\n')
                with self.assertRaises(ValueError):r.compare_tsv(a,b)
    def test_numeric_identifier_spelling_must_remain_exact(self):
        with tempfile.TemporaryDirectory() as d:
            for column in ['patient_id','sample_id','gene','configuration']:
                a=self.table(d,'a.tsv',f'{column}\tvalue\n100\t0.9\n')
                b=self.table(d,'b.tsv',f'{column}\tvalue\n100.0\t0.9\n')
                with self.assertRaises(ValueError):r.compare_tsv(a,b)
    def test_deduplicated_normalization_still_requires_independent_fresh_output(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)/'reference';workspace=Path(d)/'replay'
            original='sample_id\tlocked_original_norm_factor\tretained_library_TMM_factor\nDR144\t1.01\t0.98\n'
            self.table(root,'analysis/bulk_gate_inputs/published_retained_TMM_factors.tsv',original)
            with self.assertRaisesRegex(ValueError,'Missing freshly generated'):
                r.verify_generated_normalization(root,workspace)
            generated=self.table(workspace,'analysis/bulk_gate_normalization_factors.tsv',original)
            self.assertEqual(r.verify_generated_normalization(root,workspace)['numeric_cells_differing'],0)
            generated.write_text(original.replace('0.98','0.98001'))
            with self.assertRaises(ValueError):r.verify_generated_normalization(root,workspace)
if __name__=='__main__':unittest.main()
