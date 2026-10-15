"""Adversarial coefficient serialization and rank checks, using synthetic inputs."""
from pathlib import Path
import math
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from verify_edgeR_concordance_null import rank_unit,read_coefficients


class RefitNullTests(unittest.TestCase):
    def test_rank_ties_are_scale_invariant(self):
        a=rank_unit([1,1,2,4]);b=rank_unit([13,13,16,22])
        self.assertEqual(a,b)
        self.assertAlmostEqual(sum(v*v for v in a),1)
        self.assertEqual(a[0],a[1])

    def test_opposite_ranks_and_constant_vector(self):
        a=rank_unit([1,2,3]);b=rank_unit([3,2,1])
        self.assertAlmostEqual(math.fsum(x*y for x,y in zip(a,b)),-1)
        with self.assertRaises(ValueError):rank_unit([1,1,1])

    def read_fixture(self,text):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/"coefficients.tsv";p.write_text(text)
            return read_coefficients(p,2)

    def test_complement_columns_and_zero_coefficient(self):
        genes,columns=self.read_fixture("gene\tconfiguration_001\tconfiguration_002\na\t-1\t1\nb\t0\t-0\nc\t2\t-2\n")
        self.assertEqual(genes,["a","b","c"])
        self.assertEqual(columns,[[-1,0,2],[1,0,-2]])

    def test_bad_schema_nonfinite_or_duplicate_genes(self):
        fixtures=["gene\tconfiguration_001\na\t1\n",
                  "gene\tconfiguration_001\tconfiguration_002\na\tnan\t1\n",
                  "gene\tconfiguration_001\tconfiguration_002\na\t1\t-1\na\t2\t-2\n"]
        for text in fixtures:
            with self.assertRaises(ValueError):self.read_fixture(text)

    def test_asymmetric_complement_is_rejected(self):
        with self.assertRaises(ValueError):
            self.read_fixture("gene\tconfiguration_001\tconfiguration_002\na\t1\t-0.999\n")


if __name__=="__main__":unittest.main()
