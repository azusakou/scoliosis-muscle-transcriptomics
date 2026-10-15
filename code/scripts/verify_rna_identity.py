#!/usr/bin/env python3
"""Independently assess Somalier pair tables under predeclared RNA consistency gates.

This verifies sample coverage and the reported relatedness formula. Consistent
RNA profiles support pairing; they do not establish anatomy or clinical identity.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
from input_integrity import validate_pairs


def table(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header: {path}")
        names = [name.lstrip("#$ ") for name in reader.fieldnames]
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate header: {path}")
        rows = []
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed row: {path}")
            rows.append(dict(zip(names, row.values())))
        return rows


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assess(pairs, manifest, criteria):
    validate_pairs(manifest, expected_pairs=7)
    samples = {row["sample_id"]: row["patient_id"] for row in manifest}
    if len(samples) != len(manifest) or len(samples) != 14:
        raise ValueError("Require exactly 14 distinct manifest samples")
    patients = set(samples.values())
    if len(patients) != 7 or any(list(samples.values()).count(p) != 2 for p in patients):
        raise ValueError("Require seven two-sample patient groups")
    observed = {}
    for row in pairs:
        a, b = row["sample_a"], row["sample_b"]
        key = tuple(sorted((a, b)))
        if a == b or a not in samples or b not in samples or key in observed:
            raise ValueError(f"Invalid or duplicate comparison: {key}")
        values = {field: float(row[field]) for field in
                  ("relatedness", "ibs0", "ibs2", "n", "shared_hets", "hets_ab")}
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError(f"Nonfinite comparison: {key}")
        for field in ("ibs0", "ibs2", "n", "shared_hets", "hets_ab"):
            if values[field] < 0 or not values[field].is_integer():
                raise ValueError(f"Invalid genotype count {field}: {key}")
        if values["ibs0"] + values["ibs2"] > values["n"]:
            raise ValueError(f"IBS counts exceed shared callable sites: {key}")
        if (values["shared_hets"] > values["ibs2"] or
                2 * values["shared_hets"] > values["hets_ab"] or
                values["hets_ab"] > 2 * values["n"]):
            raise ValueError(f"Inconsistent shared-heterozygote counts: {key}")
        # Shared callable sites partition into identical genotypes (IBS2),
        # opposing homozygotes (IBS0), and het/hom discordances. Each shared
        # heterozygote contributes two to hets_ab; each discordance contributes one.
        discordances = values["n"] - values["ibs0"] - values["ibs2"]
        if values["hets_ab"] != 2 * values["shared_hets"] + discordances:
            raise ValueError(f"Inconsistent callable-site/heterozygote partition: {key}")
        expected = 2 * (values["shared_hets"] - 2 * values["ibs0"]) / max(1, values["hets_ab"])
        # Somalier 0.3.5 emits relatedness rounded to three decimal places.
        reported = float(row.get("somalier_reported_relatedness", row["relatedness"]))
        if not math.isfinite(reported) or abs(reported - expected) > 0.0005001:
            raise ValueError(f"Relatedness formula mismatch: {key}")
        if "somalier_reported_relatedness" in row and abs(values["relatedness"] - expected) > 1e-12:
            raise ValueError(f"Canonical unrounded relatedness mismatch: {key}")
        # Apply the unchanged gates to full-precision arithmetic, not a metric
        # rounded to three decimals (0.89975 must not become a passing 0.900).
        values["reported_relatedness"] = reported
        values["relatedness"] = expected
        observed[key] = values
    if set(observed) != set(itertools.combinations(sorted(samples), 2)):
        raise ValueError("Incomplete all-pairs matrix; require all 91 comparisons")
    results = []
    for patient in sorted(patients):
        a, b = sorted(s for s in samples if samples[s] == patient)
        values = observed[(a, b)]
        relatedness = values["relatedness"]
        competitors_a = [v["relatedness"] for k, v in observed.items() if a in k and b not in k]
        competitors_b = [v["relatedness"] for k, v in observed.items() if b in k and a not in k]
        reciprocal = relatedness > max(competitors_a) and relatedness > max(competitors_b)
        sufficient = values["n"] >= criteria["minimum_shared_callable_sites"]
        ibs_fraction = values["ibs0"] / values["n"] if values["n"] else None
        consistent = (sufficient and relatedness >= criteria["same_profile_relatedness_minimum"]
                      and ibs_fraction <= criteria["maximum_IBS0_fraction"]
                      and (reciprocal or not criteria["require_reciprocal_best_match"]))
        results.append({"patient_id": patient, "sample_a": a, "sample_b": b,
                        "shared_callable_sites": int(values["n"]), "relatedness": relatedness,
                        "somalier_reported_relatedness": values["reported_relatedness"],
                        "ibs0": int(values["ibs0"]), "ibs0_fraction": ibs_fraction,
                        "maximum_unpaired_relatedness_a": max(competitors_a),
                        "maximum_unpaired_relatedness_b": max(competitors_b),
                        "reciprocal_unique_best_match": reciprocal,
                        "interpretation": "consistent_RNA_profiles" if consistent else
                          ("inconclusive_coverage" if not sufficient else "does_not_meet_consistency_gates")})
    return {"structural_and_formula_checks": "PASS", "samples": len(samples),
            "all_pair_comparisons": len(observed), "expected_pairs": results,
            "scope": criteria["scope"]}


def verify_mapping_qc(rows, manifest, pairs):
    """Check full-read and callable-site denominators in completed aggregate QC."""
    by_sample = {row["sample_id"]: row for row in rows}
    if len(rows) != 14 or len(by_sample) != 14 or set(by_sample) != {r["sample_id"] for r in manifest}:
        raise ValueError("Require complete distinct 14-library mapping QC")
    site_totals = set()
    for expected in manifest:
        row = by_sample[expected["sample_id"]]
        if any(row[key] != expected[key] for key in ("patient_id", "side", "sra_run")):
            raise ValueError("Mapping-QC sample metadata mismatch")
        fields = ("validated_read_pairs", "salmon_processed_pairs", "salmon_mapped_pairs",
                  "STAR_input_pairs", "STAR_uniquely_mapped_pairs", "STAR_multimapped_pairs",
                  "somalier_called_sites", "somalier_unknown_sites", "somalier_n_het")
        counts = {key: float(row[key]) for key in fields}
        if any(not math.isfinite(v) or not v.is_integer() or v < 0 for v in counts.values()):
            raise ValueError("Invalid mapping-QC count")
        n = counts["validated_read_pairs"]
        if n <= 0 or n != counts["salmon_processed_pairs"] or n != counts["STAR_input_pairs"]:
            raise ValueError("FASTQ/Salmon/STAR input-count mismatch")
        if (counts["salmon_mapped_pairs"] > n or
                counts["STAR_uniquely_mapped_pairs"] + counts["STAR_multimapped_pairs"] > n or
                counts["somalier_n_het"] > counts["somalier_called_sites"]):
            raise ValueError("Mapping/callable-site count exceeds its denominator")
        for key, count_key, tolerance in (("salmon_percent_mapped", "salmon_mapped_pairs", 1e-8),
                                           ("STAR_uniquely_mapped_percent", "STAR_uniquely_mapped_pairs", .005001)):
            value = float(row[key])
            if not math.isfinite(value) or abs(value - 100 * counts[count_key] / n) > tolerance:
                raise ValueError("Mapping percentage disagrees with count")
        depth = float(row["somalier_genotype_depth_mean"])
        if not math.isfinite(depth) or depth < 0:
            raise ValueError("Invalid mean genotype depth")
        site_totals.add(counts["somalier_called_sites"] + counts["somalier_unknown_sites"])
    if len(site_totals) != 1 or next(iter(site_totals)) <= 0:
        raise ValueError("Fingerprint QC does not use the same site universe")
    for row in pairs:
        n = int(row["n"])
        if n > min(int(by_sample[row[s]]["somalier_called_sites"]) for s in ("sample_a", "sample_b")):
            raise ValueError("Pair callable sites exceed individual callable sites")
    return {"status": "PASS", "libraries": 14,
            "full_validated_read_pairs": sum(int(r["validated_read_pairs"]) for r in rows),
            "fingerprint_site_universe": int(next(iter(site_totals)))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("pairs", "manifest", "criteria", "output"):
        parser.add_argument("--" + argument, type=Path, required=True)
    parser.add_argument("--mapping-qc", type=Path)
    parser.add_argument("--expected-pairs", type=Path)
    args = parser.parse_args()
    criteria = json.loads(args.criteria.read_text())
    pairs, manifest = table(args.pairs), table(args.manifest)
    result = assess(pairs, manifest, criteria)
    if args.mapping_qc:
        result["full_read_mapping_QC_checks"] = verify_mapping_qc(table(args.mapping_qc), manifest, pairs)
    if args.expected_pairs:
        stored = table(args.expected_pairs)
        indexed = {row["patient_id"]: row for row in stored}
        if len(stored) != 7 or len(indexed) != 7:
            raise ValueError("Require exactly seven expected-pair aggregate rows")
        for expected in result["expected_pairs"]:
            if expected["patient_id"] not in indexed:
                raise ValueError("Expected-pair patient missing from aggregate")
            row = indexed[expected["patient_id"]]
            for key, value in expected.items():
                if isinstance(value, bool):
                    if row[key].lower() not in {str(value).lower(), str(int(value))}:
                        raise ValueError("Expected-pair boolean mismatch: " + key)
                elif isinstance(value, (int, float)):
                    number = float(row[key])
                    if not math.isfinite(number) or abs(number - value) > 1e-12:
                        raise ValueError("Expected-pair numeric mismatch: " + key)
                elif row[key] != value:
                    raise ValueError("Expected-pair interpretation/label mismatch: " + key)
        result["expected_pair_aggregate_checks"] = {"status": "PASS", "rows": 7}
        result["expected_pair_aggregate_sha256"] = sha(args.expected_pairs)
    result["prespecified_criteria"] = criteria
    result["source_sha256"] = {key: sha(getattr(args, key)) for key in ("pairs", "manifest", "criteria")}
    if args.mapping_qc:
        result["source_sha256"]["mapping_qc"] = sha(args.mapping_qc)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
