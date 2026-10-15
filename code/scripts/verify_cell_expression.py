#!/usr/bin/env python3
"""Verify actual cell-level source records without using Seurat or the R generator."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENES = {"IGKC", "ACTA1", "MYL2", "DES", "CFD", "APOD", "POSTN", "PRG4", "MFAP5", "ESR1"}


def rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def close(actual, expected, label):
    if not math.isfinite(actual) or abs(actual - float(expected)) > 2e-10 * max(1, abs(float(expected))):
        raise AssertionError(f"{label}: {actual} != {expected}")


def quantile(values, p):
    ordered = sorted(values)
    position = (len(ordered) - 1) * p  # R quantile default type 7.
    lower = math.floor(position)
    upper = math.ceil(position)
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def verify(root=ROOT):
    metadata = list(rows(root / "results/tables/annotation/annotated_cell_metadata.tsv.gz"))
    cells = {r["cell_id"]: r for r in metadata}
    assert len(cells) == len(metadata) == 14846
    assert len({r["sample_id"] for r in metadata}) == 4
    assert len({r["cell_type"] for r in metadata}) == 9
    assert all(r["passes_qc"] == "TRUE" and r["scDblFinder.class"] == "singlet" for r in metadata)
    rev = root / "analysis"
    summaries = list(rows(rev / "signature_expression_by_library_lineage_before_after.tsv"))
    summary = {(r["assay"], r["sample_id"], r["cell_type"], r["gene"]): r for r in summaries}
    assert len(summaries) == len(summary) == 11376
    groups, observed, coverage = defaultdict(list), set(), Counter()
    n = 0
    selected_totals = defaultdict(float)
    positive_records = {}
    for row in rows(rev / "signature_unexpected_gene_cell_source.tsv.gz"):
        assay, cell, gene = row["assay"], row["cell_id"], row["gene"]
        assert assay in {"RNA", "SoupX"} and gene in GENES and cell in cells
        identity = assay, cell, gene
        assert identity not in observed, f"Duplicate source record: {identity}"
        observed.add(identity)
        info = cells[cell]
        for field in ("sample_id", "patient_id", "side", "cell_type"):
            assert info[field] == row[field], f"Cell label mismatch: {cell}/{field}"
        umi = float(row["UMI"])
        assert math.isfinite(umi) and umi >= 0 and umi.is_integer()
        total = float(info["nCount_" + assay])
        assert total >= umi and total > 0
        normalized = umi / total * 10000
        ln = math.log1p(normalized)
        close(ln, row["log1p_cp10k"], f"Per-cell normalization/{identity}")
        key = assay, row["sample_id"], row["cell_type"], gene
        groups[key].append((umi, normalized, ln))
        selected_totals[(assay, cell)] += umi
        if umi > 0:
            positive_records[identity] = row
        coverage[(assay, gene)] += 1
        n += 1
    assert n == 296920
    assert set(coverage) == {(a, g) for a in ("RNA", "SoupX") for g in GENES}
    assert all(count == 14846 for count in coverage.values())
    assert all(value <= float(cells[cell]["nCount_" + assay])
               for (assay, cell), value in selected_totals.items())
    checked = 0
    for key, stored in summary.items():
        if key[-1] not in GENES:
            continue
        values = groups.get(key, [])
        assert len(values) == int(stored["n_cells"])
        if not values:
            assert int(stored["n_detected"]) == 0
            for field in ("detection_fraction", "mean_UMI_per_cell", "mean_log1p_cp10k"):
                assert stored[field] == "NA"
            continue
        umis, normalized, logvalues = zip(*values)
        detected = sum(v > 0 for v in umis)
        assert detected == int(stored["n_detected"])
        calculations = {
            "detection_fraction": detected / len(values),
            "mean_UMI_per_cell": statistics.fmean(umis),
            "mean_cp10k": statistics.fmean(normalized),
            "mean_log1p_cp10k": statistics.fmean(logvalues),
            "median_log1p_cp10k": statistics.median(logvalues),
            "p90_log1p_cp10k": quantile(logvalues, .9),
            "max_log1p_cp10k": max(logvalues),
        }
        if detected:
            calculations["mean_UMI_detected_cells"] = math.fsum(umis) / detected
        else:
            assert stored["mean_UMI_detected_cells"] == "NA"
        for field, value in calculations.items():
            close(value, stored[field], f"Aggregate/{key}/{field}")
        checked += 1
    assert checked == 700
    qc_fields = ("nCount_RNA", "nFeature_RNA", "nCount_SoupX", "nFeature_SoupX",
                 "subsets_Mito_percent", "subsets_Ribo_percent", "scDblFinder.score")
    qc_seen = set()
    for row in rows(rev / "signature_positive_cell_QC_trace.tsv.gz"):
        key = row["assay"], row["cell_id"], row["gene"]
        assert key in positive_records and key not in qc_seen
        qc_seen.add(key)
        info, source = cells[row["cell_id"]], positive_records[key]
        for field in ("UMI", "log1p_cp10k"):
            close(float(row[field]), source[field], f"Positive-cell source/{key}/{field}")
        for field in ("sample_id", "patient_id", "side", "cell_type"):
            assert row[field] == source[field]
        for field in qc_fields:
            close(float(row[field]), info[field], f"Positive-cell QC/{key}/{field}")
        for field in ("scDblFinder.class", "passes_qc_metrics", "passes_qc",
                      "seurat_clusters", "annotation_confidence"):
            assert row[field] == info[field]
    assert qc_seen == set(positive_records)
    concentration_seen = set()
    for row in rows(rev / "signature_positive_cell_QC_concentration.tsv"):
        key = row["assay"], row["sample_id"], row["cell_type"], row["gene"]
        assert key in groups and key not in concentration_seen
        concentration_seen.add(key)
        values = groups[key]
        positive = [v[0] for v in values if v[0] > 0]
        infos = [cells[cell] for (assay, cell, gene), source in positive_records.items()
                 if (assay, source["sample_id"], source["cell_type"], gene) == key]
        assert len(values) == int(row["n_cells"]) and len(positive) == int(row["n_detected"])
        close(sum(positive), row["total_gene_UMI"], f"UMI total/{key}")
        close(max(positive, default=0), row["largest_single_cell_UMI"], f"Largest UMI/{key}")
        if positive:
            calculations = {
                "largest_cell_fraction_of_gene_UMIs": max(positive) / sum(positive),
                "positive_cells_median_nCount_RNA": statistics.median(float(i["nCount_RNA"]) for i in infos),
                "positive_cells_max_mito_percent": max(float(i["subsets_Mito_percent"]) for i in infos),
                "positive_cells_max_original_doublet_score": max(float(i["scDblFinder.score"]) for i in infos)}
            for field, value in calculations.items():
                close(value, row[field], f"Positive-cell concentration/{key}/{field}")
            assert row["all_positive_cells_original_QC_singlets"] == "TRUE"
        else:
            assert row["largest_cell_fraction_of_gene_UMIs"] == "NA"
    assert concentration_seen == set(groups)
    return {"status": "PASS", "cells": len(cells), "source_records": n,
            "independently_recomputed_nonempty_groups": checked,
            "independently_verified_positive_QC_records": len(qc_seen),
            "independently_recomputed_QC_concentration_groups": len(concentration_seen),
            "scope": "Ten nominated genes; both assays; labels, counts, per-cell normalization and distribution summaries"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = verify()
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
