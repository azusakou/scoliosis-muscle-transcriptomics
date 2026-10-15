import copy
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import verify_remaining_evidence as v
import verify_test_family_inventory as inventory


class SourceLinkTests(unittest.TestCase):
    def test_proxy_patient_identity_and_expression_are_source_bound(self):
        original = v.rows('bulk_followup_ADIPOQ_proxy_patient_effects.tsv')
        v.verify_proxy_input_links(original)
        for field,value in [('patient_id','wrong-patient'),('ADIPOQ_difference','100'),('cohort','wrong-cohort')]:
            changed = copy.deepcopy(original);changed[0][field] = value
            with self.assertRaises((ValueError,AssertionError)):v.verify_proxy_input_links(changed)

    def test_clinical_covariate_score_and_identity_join_reject_corruption(self):
        original = v.rows('bulk_followup_available_clinical_patient_effects.tsv')
        v.verify_clinical_input_links(original)
        for field,value in [('age_years','99'),('score_difference','100'),('analysis_patient_id','wrong-patient')]:
            changed = copy.deepcopy(original);changed[0][field] = value
            with self.assertRaises((ValueError,AssertionError)):v.verify_clinical_input_links(changed)

    def test_composition_method_patient_lineage_and_values_are_source_bound(self):
        original = v.rows('bulk_followup_composition_patient_differences.tsv')
        v.verify_composition_input_links(original)
        for field,value in [('patient_id','wrong-patient'),('method','wrong-method'),('B/plasma cell','1')]:
            changed = copy.deepcopy(original);changed[0][field] = value
            with self.assertRaises((ValueError,AssertionError)):v.verify_composition_input_links(changed)

    def test_inventory_rejects_a_missing_signature_test_block(self):
        original = v.rows('bulk_test_family_inventory.tsv')
        inventory.validate(v.ROOT, original)
        changed = [r for r in original if r['source'] != 'analysis/signature_bulk_sensitivity_tests.tsv']
        with self.assertRaisesRegex(ValueError,'Omitted P-value artifacts'):inventory.validate(v.ROOT, changed)


class InventoryCoverageTests(unittest.TestCase):
    def test_raw_probability_column_suffixes_are_recognized(self):
        columns=['PValue.AIS','PValue.CS','PValue.muscle_GSE254300','published_pvalue','meta_p','exact_signflip_p','p_val']
        self.assertEqual(inventory.p_columns(columns), columns)
        self.assertEqual(inventory.p_columns(['FDR.AIS','adj.P.Val','n_permutations','minimum_attainable_two_sided_P','p_reproduction_error','reported_P_columns']), [])

    def test_dot_qualified_gene_comparison_artifacts_cannot_be_omitted(self):
        records=v.rows('bulk_test_family_inventory.tsv')
        sources=['results/tables/adjacent_tissue_gse324067/muscle_osteoblast_effect_comparison.tsv.gz',
                 'results/tables/disease_control_gse254300_cs/GSE254300_AIS_CS_gene_effect_comparison.tsv.gz',
                 'results/tables/published_summary_gse311180/GSE254300_GSE311180_effect_comparison.tsv.gz']
        for source in sources:
            with self.subTest(source=source):
                with self.assertRaisesRegex(ValueError,'Omitted P-value artifacts'):
                    inventory.validate(v.ROOT,[r for r in records if r['source'] != source])

    def test_both_figure_source_directories_must_be_covered(self):
        records=v.rows('bulk_test_family_inventory.tsv')
        for source in ['results/figures/source_data/Figure_S4_model_sensitivity.tsv',
                       'analysis/tables/figure_validation_t_confidence_intervals.tsv']:
            with self.subTest(source=source):
                with self.assertRaisesRegex(ValueError,'Omitted P-value artifacts'):
                    inventory.validate(v.ROOT,[r for r in records if r['source'] != source])

if __name__ == '__main__':unittest.main()

