import csv, math, sys, tempfile, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import verify_results as verify

class TemporaryCase(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
    def tearDown(self):
        self.temp.cleanup()

class StatisticalTests(TemporaryCase):
    def test_t_distribution_known_values(self):
        for value in [0, .2, 1, 4, 20]:
            self.assertAlmostEqual(verify.t_two_sided(value, 1), 1 - 2 * math.atan(abs(value)) / math.pi, places=12)
        self.assertAlmostEqual(verify.t_critical(11), 2.200985160091639, places=11)
    def test_complete_sign_flip_and_wilcoxon_tie_guard(self):
        self.assertEqual(verify.exact_side([1, 2, 3]), .25)
        self.assertEqual(verify.exact_side([0, 1]), 1)
        self.assertEqual(verify.exact_wilcoxon([1, 2, 3]), (6, .25))
        with self.assertRaises(AssertionError):
            verify.exact_wilcoxon([1, 1, -1])
    def test_wrong_bh_output_is_rejected(self):
        file = self.root / 'results/tables/invalid.tsv'; file.parent.mkdir(parents=True)
        file.write_text('PValue\tFDR\n0.01\t0.9\n0.5\t0.5\n')
        original_root = verify.ROOT
        try:
            verify.ROOT = self.root
            with self.assertRaises(AssertionError):
                verify.schema_and_fdr()
        finally:
            verify.ROOT = original_root
    def test_bh_ties_and_empty(self):
        self.assertEqual(verify.bh([]), [])
        self.assertEqual(verify.bh([.01, .01, .04, .9]), [.02, .02, .05333333333333334, .9])
    def test_broken_tsv_is_rejected(self):
        file = self.root / 'bad.tsv'; file.write_text('a\tb\n1\n')
        with self.assertRaises(AssertionError):
            verify.read_table(file)
        file.write_text('a\tb\n"unclosed\t2\n')
        with self.assertRaises(csv.Error):
            verify.read_table(file)
        file.write_text('a\ta\n1\t2\n')
        with self.assertRaises(AssertionError):
            verify.read_table(file)
