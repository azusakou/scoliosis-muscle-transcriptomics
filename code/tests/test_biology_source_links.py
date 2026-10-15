import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_biology_cells import source_value, aggregate


class BiologicalSourceTests(unittest.TestCase):
    def setUp(self):
        self.info = dict(sample_id="D1__concave", patient_id="D1", side="concave",
                         cell_type="FAP", nCount_RNA="100", nCount_SoupX="80")
        self.row = dict(assay="SoupX", cell_id="c1", gene="ANGPT1", UMI="2",
                        all_gene_total_UMI="80", cp10k="250", log1p_cp10k=str(math.log1p(250)),
                        **{f: self.info[f] for f in ("sample_id", "patient_id", "side", "cell_type")})

    def test_wrong_patient_or_lineage_cannot_pass_arithmetic(self):
        source_value(self.row, {"c1": self.info})
        for field in ("patient_id", "cell_type", "side"):
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, "label mismatch"):
                source_value({**self.row, field: "wrong"}, {"c1": self.info})

    def test_wrong_library_total_rejected_even_if_normalized_values_match(self):
        changed = {**self.row, "all_gene_total_UMI": "100", "cp10k": "200",
                   "log1p_cp10k": str(math.log1p(200))}
        with self.assertRaisesRegex(ValueError, "source library total"):
            source_value(changed, {"c1": self.info})

    def test_empty_group_distinct_from_detected_zero_expression(self):
        empty, zero = aggregate([]), aggregate([(0., 0., 0.)])
        self.assertIsNone(empty["detection_fraction"])
        self.assertEqual(zero["detection_fraction"], 0.)
        self.assertIsNone(zero["fraction_UMI_in_top1_cell"])

    def test_sparse_concentration_retains_zero_cells_in_denominator(self):
        z = aggregate([(10., 100., math.log1p(100))] + [(0., 0., 0.)]*99)
        self.assertEqual(z["detection_fraction"], .01)
        self.assertEqual(z["fraction_UMI_in_top1percent_cells"], 1.)
        self.assertEqual(z["mean_UMI_per_cell"], .1)


if __name__ == "__main__":
    unittest.main()
