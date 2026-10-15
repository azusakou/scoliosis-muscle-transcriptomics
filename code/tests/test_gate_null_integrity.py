import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_gate_concordance import exact_grid
from verify_test_family_inventory import validate


class GateNullIntegrityTests(unittest.TestCase):
    def test_duplicate_configuration_is_rejected(self):
        rows = [dict(configuration=str(i+1), **{"p"+str(j): str(1 if (i >> j)&1 else -1) for j in range(5)}) for i in range(32)]
        exact_grid(rows, 5)
        rows[5] = dict(rows[4])
        with self.assertRaisesRegex(ValueError, "configuration"):
            exact_grid(rows, 5)

    def test_missing_configuration_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "grid size"):
            exact_grid([dict(configuration="1", p0="-1", p1="-1")], 2)

    def test_inventory_must_cover_compressed_gene_refits(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); rev = root/"analysis"; rev.mkdir(parents=True)
            (root/"results/tables").mkdir(parents=True)
            path = rev/"bulk_gate_genes.tsv.gz"
            with gzip.open(path, "wt") as f: f.write("gene\tPValue\tFDR\nA\t0.01\t0.02\n")
            with self.assertRaisesRegex(ValueError, "Omitted P-value artifacts"):
                validate(root, [])
            result = validate(root, [dict(family="fixed_refit",source="analysis/bulk_gate_genes.tsv.gz",
                                           role="gene_refit",correction="BH_across_fixed_gene_universe")])
            self.assertEqual(result["P_value_artifacts_covered"], 1)


if __name__ == "__main__": unittest.main()
