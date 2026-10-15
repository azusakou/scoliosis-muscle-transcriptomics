"""Adversarial checks for interpretation of the complete RNA pair matrix."""
import copy
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_rna_identity import assess,verify_mapping_qc


class RnaIdentityTests(unittest.TestCase):
    def setUp(self):
        self.manifest = [{"sample_id": f"P{p}_{s}", "patient_id": f"P{p}",
                          "side": "concave" if s=="a" else "convex"}
                         for p in range(7) for s in ("a", "b")]
        self.criteria = {"minimum_shared_callable_sites": 1000,
                         "same_profile_relatedness_minimum": .9,
                         "maximum_IBS0_fraction": .01,
                         "require_reciprocal_best_match": True,
                         "scope": "RNA pairing only"}
        self.pairs = []
        for a, b in itertools.combinations(self.manifest, 2):
            same = a["patient_id"] == b["patient_id"]
            self.pairs.append({"sample_a": a["sample_id"], "sample_b": b["sample_id"],
                               "n": "2000", "ibs0": "0", "ibs2": "2000" if same else "1100",
                               "hets_ab": "1000", "shared_hets": "500" if same else "50",
                               "relatedness": "1.000" if same else "0.100"})

    def test_complete_consistent_profiles(self):
        result = assess(self.pairs, self.manifest, self.criteria)
        self.assertEqual(result["all_pair_comparisons"], 91)
        self.assertTrue(all(r["interpretation"] == "consistent_RNA_profiles"
                            for r in result["expected_pairs"]))

    def test_missing_or_duplicated_comparisons_rejected(self):
        for pairs in (self.pairs[:-1], self.pairs + [self.pairs[0]]):
            with self.assertRaises(ValueError):
                assess(pairs, self.manifest, self.criteria)

    def test_invented_statistic_rejected(self):
        self.pairs[0]["relatedness"] = "0.998"
        with self.assertRaises(ValueError):
            assess(self.pairs, self.manifest, self.criteria)

    def test_ties_do_not_establish_unique_matching(self):
        self.pairs[1]["relatedness"] = "1.000"
        self.pairs[1]["shared_hets"] = "500"
        self.pairs[1]["ibs2"] = "2000"
        result = assess(self.pairs, self.manifest, self.criteria)
        self.assertEqual(result["expected_pairs"][0]["interpretation"],
                         "does_not_meet_consistency_gates")

    def test_low_coverage_is_inconclusive(self):
        self.pairs[0]["n"] = "999"
        self.pairs[0]["ibs2"] = "999"
        result = assess(self.pairs, self.manifest, self.criteria)
        self.assertEqual(result["expected_pairs"][0]["interpretation"], "inconclusive_coverage")

    def test_rounded_metric_does_not_cross_predeclared_gate(self):
        self.pairs[0].update(relatedness="0.900",shared_hets="3599",hets_ab="8000",
                             n="16000",ibs2="15198")
        result=assess(self.pairs,self.manifest,self.criteria)
        self.assertAlmostEqual(result["expected_pairs"][0]["relatedness"],.89975)
        self.assertEqual(result["expected_pairs"][0]["interpretation"],"does_not_meet_consistency_gates")

    def test_fraction_and_malformed_counts(self):
        for field, value in (("n", "nan"), ("ibs0", "-1"), ("ibs2", "2001"), ("hets_ab", "2.5")):
            pairs = copy.deepcopy(self.pairs)
            pairs[0][field] = value
            with self.assertRaises(ValueError):
                assess(pairs, self.manifest, self.criteria)

    def test_canonical_preserves_native_display_and_full_precision(self):
        self.pairs[0].update(relatedness=str(2*482/999),somalier_reported_relatedness="0.965",
                             shared_hets="482",hets_ab="999",ibs2="1965")
        result=assess(self.pairs,self.manifest,self.criteria)
        self.assertEqual(result["expected_pairs"][0]["somalier_reported_relatedness"],.965)
        self.assertEqual(result["expected_pairs"][0]["relatedness"],2*482/999)

    def test_corrupt_canonical_full_precision_is_rejected(self):
        self.pairs[0].update(relatedness="0.9999",somalier_reported_relatedness="1.000")
        with self.assertRaises(ValueError):assess(self.pairs,self.manifest,self.criteria)

    def test_impossible_heterozygote_counts_are_rejected(self):
        for changes in ({"ibs2":"499"},{"hets_ab":"999"},{"hets_ab":"4001"}):
            pairs=copy.deepcopy(self.pairs);pairs[0].update(changes)
            with self.assertRaises(ValueError):assess(pairs,self.manifest,self.criteria)

    def test_formula_consistent_but_impossible_site_partition_is_rejected(self):
        # The reported relatedness can agree with its formula while the
        # underlying six counts cannot describe any genotype-pair table.
        for changes in ({"ibs2": "1999"}, {"n": "2001"},
                        {"ibs0": "1", "relatedness": "0.996"}):
            pairs = copy.deepcopy(self.pairs)
            pairs[0].update(changes)
            with self.assertRaisesRegex(ValueError, "partition|IBS counts"):
                assess(pairs, self.manifest, self.criteria)

    def test_two_same_side_libraries_are_not_a_valid_pair(self):
        self.manifest[0]["side"]="convex"
        with self.assertRaises(ValueError):assess(self.pairs,self.manifest,self.criteria)

    def mapping_fixture(self):
        rows=[]
        for sample in self.manifest:
            sample["sra_run"]="synthetic_"+sample["sample_id"]
            rows.append({**sample,"validated_read_pairs":"100","salmon_processed_pairs":"100",
                         "salmon_mapped_pairs":"95","salmon_percent_mapped":"95",
                         "STAR_input_pairs":"100","STAR_uniquely_mapped_pairs":"80",
                         "STAR_multimapped_pairs":"20","STAR_uniquely_mapped_percent":"80",
                         "somalier_called_sites":"2000","somalier_unknown_sites":"500",
                         "somalier_n_het":"1000","somalier_genotype_depth_mean":"20"})
        return rows

    def test_full_read_mapping_denominators(self):
        result=verify_mapping_qc(self.mapping_fixture(),self.manifest,self.pairs)
        self.assertEqual(result["full_validated_read_pairs"],1400)
        self.assertEqual(result["fingerprint_site_universe"],2500)

    def test_incomplete_reads_and_inconsistent_QC_are_rejected(self):
        rows=self.mapping_fixture()
        for key,value in (("STAR_input_pairs","99"),("salmon_percent_mapped","94.9"),
                          ("STAR_uniquely_mapped_percent","80.006"),
                          ("somalier_called_sites","1999"),("somalier_genotype_depth_mean","nan")):
            altered=copy.deepcopy(rows);altered[0][key]=value
            with self.assertRaises(ValueError):verify_mapping_qc(altered,self.manifest,self.pairs)


if __name__ == "__main__":
    unittest.main()
