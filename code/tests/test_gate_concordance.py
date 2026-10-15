import copy
import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_gate_concordance import (
    REQUIRED_GATE_INPUTS, verify_input_bindings, verify_normalization_join,
)


class FixedNormalizationTests(unittest.TestCase):
    def setUp(self):
        self.original = [
            dict(sample_id="P1c", patient_id="P1", side="concave"),
            dict(sample_id="P1v", patient_id="P1", side="convex"),
            dict(sample_id="P2c", patient_id="P2", side="concave"),
            dict(sample_id="P2v", patient_id="P2", side="convex"),
        ]
        self.qc = [dict(row, normalization_factor=str(value)) for row, value in
                   zip(self.original, (1.01, .97, 1.12, 1.05))]
        self.metadata = copy.deepcopy(self.original[:2])
        self.published = [dict(row, locked_original_norm_factor=self.qc[i]["normalization_factor"],
                               retained_library_TMM_factor=str((1.02, .98)[i]),
                               retained_tested_library_size=str((1000., 1100.)[i]))
                          for i, row in enumerate(self.metadata)]
        self.reported = copy.deepcopy(self.published)

    def check(self, qc=None, published=None, reported=None, metadata=None, original=None):
        return verify_normalization_join(
            self.qc if qc is None else qc, self.published if published is None else published,
            self.reported if reported is None else reported,
            self.metadata if metadata is None else metadata,
            self.original if original is None else original)

    def test_join_uses_sample_identity_rather_than_row_position(self):
        result = self.check(qc=list(reversed(self.qc)), published=list(reversed(self.published)))
        self.assertEqual(result["locked_factors_joined"], 2)
        self.assertEqual(result["published_retained_factors_joined"], 2)
        self.assertEqual(result["original_QC_sample_ids_verified"], 4)

    def test_locked_factor_cannot_change_even_if_published_copy_matches(self):
        p, r = copy.deepcopy(self.published), copy.deepcopy(self.reported)
        p[0]["locked_original_norm_factor"] = r[0]["locked_original_norm_factor"] = "1.010001"
        with self.assertRaisesRegex(ValueError, "locked factor QC join"):
            self.check(published=p, reported=r)

    def test_native_retained_factor_drift_cannot_replace_published_parameter(self):
        r = copy.deepcopy(self.reported)
        r[0]["retained_library_TMM_factor"] = str(1.02 + 5.618824250319676e-7)
        with self.assertRaisesRegex(ValueError, "published retained factor join"):
            self.check(reported=r)

    def test_matching_factors_do_not_override_wrong_patient_or_side(self):
        for field, value in (("patient_id", "P2"), ("side", "convex")):
            with self.subTest(field=field):
                p = copy.deepcopy(self.published); p[0][field] = value
                with self.assertRaisesRegex(ValueError, "identity mismatch"):
                    self.check(published=p)

    def test_duplicate_sample_ids_rejected_in_every_join_source(self):
        for name in ("qc", "published", "reported", "metadata", "original"):
            with self.subTest(source=name):
                rows = copy.deepcopy(getattr(self, name)); rows.append(copy.deepcopy(rows[0]))
                with self.assertRaisesRegex(ValueError, "duplicate normalization sample ID"):
                    self.check(**{name: rows})

    def test_missing_or_unexpected_parameter_samples_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "sample coverage"):
            self.check(published=self.published[:1])
        p = copy.deepcopy(self.published); p[0]["sample_id"] = "unexpected"
        with self.assertRaisesRegex(ValueError, "sample coverage"):
            self.check(published=p)

    def test_nonpositive_and_nonfinite_values_rejected_in_all_parameter_sources(self):
        fields = (("qc", "normalization_factor"),
                  ("published", "locked_original_norm_factor"),
                  ("published", "retained_library_TMM_factor"),
                  ("reported", "locked_original_norm_factor"),
                  ("reported", "retained_library_TMM_factor"),
                  ("reported", "retained_tested_library_size"))
        for source, field in fields:
            for bad in ("0", "-1", "nan", "inf", "NA"):
                with self.subTest(source=source, field=field, value=bad):
                    rows = copy.deepcopy(getattr(self, source)); rows[0][field] = bad
                    with self.assertRaisesRegex(ValueError, "Nonpositive/nonfinite"):
                        self.check(**{source: rows})

    def test_serialization_epsilon_uses_existing_arithmetic_tolerance(self):
        r = copy.deepcopy(self.reported)
        r[0]["locked_original_norm_factor"] = str(1.01 + 1e-15)
        r[0]["retained_library_TMM_factor"] = str(1.02 - 1e-15)
        self.check(reported=r)


class FactorInputBindingTests(unittest.TestCase):
    def fixtures(self, root):
        rows = []
        for name in sorted(REQUIRED_GATE_INPUTS):
            data = ("immutable input " + name).encode()
            (root/name).write_bytes(data)
            rows.append(dict(artifact=name, SHA256=hashlib.sha256(data).hexdigest()))
        return rows

    def manifest(self, root, rows):
        with (root/"source_bindings.tsv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=("artifact", "SHA256"), delimiter="\t")
            writer.writeheader(); writer.writerows(rows)

    def test_published_factor_cannot_be_omitted_from_source_bindings(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); rows = self.fixtures(root); self.manifest(root, rows)
            self.assertEqual(len(verify_input_bindings(root)), len(REQUIRED_GATE_INPUTS))
            self.manifest(root, [r for r in rows if r["artifact"] != "published_retained_TMM_factors.tsv"])
            with self.assertRaisesRegex(ValueError, "Unbound required gate inputs"):
                verify_input_bindings(root)

    def test_changed_published_factor_bytes_fail_sha_binding(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); self.manifest(root, self.fixtures(root))
            p = root/"published_retained_TMM_factors.tsv"; p.write_bytes(p.read_bytes()+b"changed")
            with self.assertRaisesRegex(ValueError, "hash mismatch: published_retained_TMM_factors"):
                verify_input_bindings(root)


if __name__ == "__main__":
    unittest.main()
