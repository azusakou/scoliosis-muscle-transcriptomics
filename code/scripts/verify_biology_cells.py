#!/usr/bin/env python3
"""Independently reconstruct the closed biological target audit from real cells.

Metadata/total-UMI joins, complete coverage and all unfavorable distributions
are checked. This does not certify anatomical side or contamination-free cells.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGETS = {"ANGPT1", "LEP", "ADIPOQ", "PDGFRA", "COL1A1", "COL1A2", "COL3A1",
           "DCN", "LUM", "PI16", "DPP4", "CFD", "APOD", "POSTN", "IGKC", "ACTA1", "MYL2", "DES"}


def rows(path):
    with (gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else path.open(newline="")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed biological-audit table")
            yield row


def close(actual, stored, label):
    if actual is None:
        if stored != "NA": raise ValueError("Expected undefined value: " + label)
        return
    expected = float(stored)
    if not math.isfinite(expected) or not math.isclose(actual, expected, rel_tol=2e-10, abs_tol=2e-12):
        raise ValueError("Biological-audit value mismatch: " + label)


def source_value(row, metadata):
    if row["assay"] not in {"RNA", "SoupX"} or row["gene"] not in TARGETS:
        raise ValueError("Unexpected assay or target")
    if row["cell_id"] not in metadata: raise ValueError("Unknown cell identity")
    info = metadata[row["cell_id"]]
    for field in ("sample_id", "patient_id", "side", "cell_type"):
        if row[field] != info[field]: raise ValueError("Cell source label mismatch: " + field)
    umi, total = float(row["UMI"]), float(row["all_gene_total_UMI"])
    if not (math.isfinite(umi) and umi >= 0 and umi.is_integer() and total > 0 and total >= umi):
        raise ValueError("Invalid count/total")
    close(total, info["nCount_" + row["assay"]], "source library total")
    normalized = umi / total * 10000
    ln = math.log1p(normalized)
    close(normalized, row["cp10k"], "cell cp10k")
    close(ln, row["log1p_cp10k"], "cell log1p cp10k")
    return umi, normalized, ln


def aggregate(values):
    n = len(values)
    if not n:
        return {"n_cells": 0, "n_detected": 0, "total_UMI": 0,
                **{f: None for f in ("detection_fraction", "mean_UMI_per_cell", "mean_cp10k",
                    "mean_log1p_cp10k", "median_log1p_cp10k", "fraction_UMI_in_top1_cell",
                    "fraction_UMI_in_top5_cells", "fraction_UMI_in_top1percent_cells")}}
    umis, normalized, logs = zip(*values)
    total = math.fsum(umis); ordered = sorted(umis, reverse=True)
    return dict(n_cells=n, n_detected=sum(v > 0 for v in umis), total_UMI=total,
                detection_fraction=sum(v > 0 for v in umis) / n,
                mean_UMI_per_cell=total / n, mean_cp10k=statistics.fmean(normalized),
                mean_log1p_cp10k=statistics.fmean(logs), median_log1p_cp10k=statistics.median(logs),
                fraction_UMI_in_top1_cell=ordered[0] / total if total else None,
                fraction_UMI_in_top5_cells=math.fsum(ordered[:5]) / total if total else None,
                fraction_UMI_in_top1percent_cells=math.fsum(ordered[:max(1, math.ceil(n*.01))]) / total if total else None)


def verify(root=ROOT):
    rev = root / "analysis"
    metadata_rows = list(rows(root / "results/tables/annotation/annotated_cell_metadata.tsv.gz"))
    metadata = {r["cell_id"]: r for r in metadata_rows}
    if len(metadata) != len(metadata_rows) or len(metadata) != 14846: raise ValueError("Cell metadata coverage")
    panel = list(rows(rev / "signature_biology_target_panel_manifest.tsv"))
    if {r["gene"] for r in panel} != TARGETS or len(panel) != len(TARGETS): raise ValueError("Closed target panel changed")
    groups, seen, shared_source = defaultdict(list), set(), {}
    prior_targets = {"IGKC", "ACTA1", "MYL2", "DES", "CFD", "APOD", "POSTN"}
    for row in rows(rev / "signature_biology_all_target_cell_source.tsv.gz"):
        identity = row["assay"], row["cell_id"], row["gene"]
        if identity in seen: raise ValueError("Duplicate cell/assay/gene")
        seen.add(identity)
        key = row["assay"], row["sample_id"], row["cell_type"], row["gene"]
        value = source_value(row, metadata)
        groups[key].append(value)
        if row["gene"] in prior_targets: shared_source[identity] = value[0]
    if len(seen) != len(metadata) * len(TARGETS) * 2: raise ValueError("Incomplete real-cell coverage")
    prior_seen = set()
    for row in rows(rev / "signature_unexpected_gene_cell_source.tsv.gz"):
        if row["gene"] not in prior_targets: continue
        key = row["assay"], row["cell_id"], row["gene"]
        if key in prior_seen or key not in shared_source: raise ValueError("Prior source linkage")
        prior_seen.add(key); close(shared_source[key], row["UMI"], "prior real-cell source count")
    if prior_seen != set(shared_source): raise ValueError("Incomplete prior source linkage")
    summaries = list(rows(rev / "signature_biology_target_gene_expression_all_libraries_lineages.tsv"))
    stored = {(r["assay"], r["sample_id"], r["cell_type"], r["gene"]): r for r in summaries}
    libraries = {r["sample_id"] for r in metadata_rows}; lineages = {r["cell_type"] for r in metadata_rows}
    expected = {(a, lib, ct, g) for a in ("RNA", "SoupX") for lib in libraries for ct in lineages for g in TARGETS}
    if set(stored) != expected or len(stored) != len(summaries): raise ValueError("Summary group coverage")
    computed = {}
    for key, r in stored.items():
        if r["feature_available"] != "TRUE": raise ValueError("Target availability changed")
        values = aggregate(groups.get(key, [])); computed[key] = values
        for field, value in values.items(): close(value, r[field], str(key) + "/" + field)
        if r["sample_id"] != r["patient_id"] + "__" + r["side"]: raise ValueError("Summary source map")
    pairs = list(rows(rev / "signature_biology_target_gene_descriptive_within_donor_differences.tsv"))
    pair_seen = set()
    for r in pairs:
        key = r["assay"], r["cell_type"], r["gene"], r["patient_id"]
        if key in pair_seen: raise ValueError("Duplicate donor difference")
        pair_seen.add(key)
        cc = computed[(key[0], key[3]+"__concave", key[1], key[2])]
        cv = computed[(key[0], key[3]+"__convex", key[1], key[2])]
        for side, z in (("concave", cc), ("convex", cv)):
            close(z["n_cells"], r["n_"+side], "paired group size")
            close(z["mean_UMI_per_cell"], r[side+"_mean_UMI"], "paired source mean")
        for source, dest in (("mean_UMI_per_cell", "mean_UMI_difference"),
                             ("mean_log1p_cp10k", "mean_log1p_cp10k_difference"),
                             ("detection_fraction", "detection_fraction_difference")):
            value = None if cc[source] is None or cv[source] is None else cc[source]-cv[source]
            close(value, r[dest], "donor descriptive difference")
    if len(pair_seen) != 2*2*len(TARGETS)*len(lineages): raise ValueError("Incomplete donor contrasts")
    shares = list(rows(rev / "signature_biology_retained_cell_expression_lineage_provenance.tsv"))
    share_seen = set()
    for r in shares:
        key = r["assay"], r["sample_id"], r["cell_type"], r["gene"]
        if key in share_seen or key not in computed: raise ValueError("Lineage share coverage")
        share_seen.add(key); z = computed[key]
        for field in ("n_cells", "n_detected", "total_UMI"): close(z[field], r[field], "lineage source")
        for source, dest in (("total_UMI", "fraction_gene_UMI_in_lineage"),
                             ("n_detected", "fraction_detected_cells_in_lineage")):
            denominator = sum(computed[(key[0], key[1], ct, key[3])][source] for ct in lineages)
            close(z[source]/denominator if denominator else None, r[dest], "lineage share")
    if share_seen != expected: raise ValueError("Incomplete lineage shares")
    return dict(status="PASS", real_cells=len(metadata), fixed_targets=len(TARGETS),
                source_cell_gene_assay_records=len(seen), reconstructed_summary_groups=len(stored),
                independently_linked_prior_source_records=len(prior_seen),
                donor_descriptive_contrasts=len(pair_seen), lineage_share_groups=len(shares),
                scope="Real exported UMI records and retained metadata; no new biological replicates or cell-level significance")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--report", type=Path); args = p.parse_args(); result = verify(args.root)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
