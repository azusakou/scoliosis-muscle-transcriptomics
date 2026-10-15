#!/usr/bin/env python3
"""Recompute all 4096 edgeR-coefficient concordance values from exported refits."""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import math
from pathlib import Path
from project_paths import resolve_analysis_input

ROOT = Path(__file__).resolve().parents[1]


def table(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed TSV: {path}")
            yield row


def close(actual, expected, label, tolerance=2e-10):
    a, b = float(actual), float(expected)
    if not math.isfinite(a) or not math.isfinite(b) or abs(a - b) > tolerance * max(1, abs(b)):
        raise ValueError(f"Numerical mismatch: {label}: {a} != {b}")


def rank_unit(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        rank = (start + end + 1) / 2
        for i in order[start:end]: ranks[i] = rank
        start = end
    center = (len(values) + 1) / 2
    shifted = [r - center for r in ranks]
    norm = math.sqrt(math.fsum(r * r for r in shifted))
    if norm == 0: raise ValueError("Constant coefficient vector")
    return [r / norm for r in shifted]


def read_coefficients(path, count):
    keys = [f"configuration_{i:03d}" for i in range(1, count + 1)]
    genes, seen, columns = [], set(), [[] for _ in keys]
    for row in table(path):
        if set(row) != {"gene", *keys} or not row["gene"] or row["gene"] in seen:
            raise ValueError(f"Invalid coefficient schema/gene: {path}")
        genes.append(row["gene"]); seen.add(row["gene"])
        for values, key in zip(columns, keys):
            value = float(row[key])
            if not math.isfinite(value): raise ValueError("Nonfinite refit coefficient")
            values.append(value)
    if not genes: raise ValueError("Empty coefficient matrix")
    for i in range(count // 2):
        for a, b in zip(columns[i], columns[-1-i]):
            close(a, -b, "complement coefficient symmetry", tolerance=1e-13)
    return genes, columns


def verify(root=ROOT):
    rev = root / "analysis"
    ga, a = read_coefficients(resolve_analysis_input("analysis/bulk_edgeR_null_logFC_GSE254300.tsv.gz", root), 32)
    gb, b = read_coefficients(rev / "bulk_edgeR_null_logFC_GSE311180.tsv.gz", 128)
    replication = list(table(root / "results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz"))
    if ga != gb or ga != [r["gene"] for r in replication]:
        raise ValueError("Refit/replication common-gene ordering mismatch")
    for j, row in enumerate(replication):
        close(a[-1][j], row["logFC.GSE254300"], "observed GSE254300 refit", 1e-7)
        close(b[-1][j], row["logFC.GSE311180_requant"], "observed GSE311180 refit", 1e-7)
    for cohort, count, patients in [("GSE254300",32,5),("GSE311180",128,7)]:
        grid = list(table(resolve_analysis_input(f"analysis/bulk_edgeR_null_grid_{cohort}.tsv", root)))
        if len(grid) != count or len(grid[0]) != patients + 1:
            raise ValueError("Invalid side-configuration grid")
        names = [key for key in grid[0] if key != "configuration"]
        for i, row in enumerate(grid):
            if int(row["configuration"]) != i + 1:
                raise ValueError("Misordered configuration")
            expected = tuple(1 if (i >> j) & 1 else -1 for j in range(patients))
            if tuple(int(row[key]) for key in names) != expected:
                raise ValueError("Patient-side sign grid mismatch")
    ra, rb = [rank_unit(v) for v in a], [rank_unit(v) for v in b]
    signs_a = [[(v > 0) - (v < 0) for v in col] for col in a]
    signs_b = [[(v > 0) - (v < 0) for v in col] for col in b]
    observed_rho = math.fsum(x*y for x,y in zip(ra[-1],rb[-1]))
    observed_agree = sum(x==y for x,y in zip(signs_a[-1],signs_b[-1])) / len(ga)
    null = list(table(rev / "bulk_gene_concordance_edgeR_permutation_null.tsv"))
    seen, rhos, agreements = set(), [], []
    for row in null:
        i,j = int(row["GSE254300_configuration"])-1,int(row["GSE311180_configuration"])-1
        if not (0 <= i < 32 and 0 <= j < 128) or (i,j) in seen:
            raise ValueError("Invalid/duplicate joint configuration")
        seen.add((i,j))
        rho = math.fsum(x*y for x,y in zip(ra[i],rb[j]))
        agree = sum(x==y for x,y in zip(signs_a[i],signs_b[j])) / len(ga)
        close(rho,row["spearman_rho"],"complete original-coefficient rho null")
        close(agree,row["sign_concordance"],"complete original-coefficient sign null")
        rhos.append(rho); agreements.append(agree)
    if len(seen) != 4096: raise ValueError("Incomplete null: require all 4096 configurations")
    summaries = list(table(rev / "bulk_gene_concordance_edgeR_permutation_summary.tsv"))
    if len(summaries) != 1: raise ValueError("Ambiguous null summary")
    summary = summaries[0]
    rho_p = sum(abs(v)>=abs(observed_rho)-1e-12 for v in rhos)/4096
    sign_p = sum(abs(v-.5)>=abs(observed_agree-.5)-1e-12 for v in agreements)/4096
    for field,value in [("observed_spearman_rho",observed_rho),("observed_sign_concordance",observed_agree),
                        ("exact_two_sided_rho_p",rho_p),("exact_two_sided_sign_concordance_p",sign_p),
                        ("n_common_genes",len(ga)),("n_joint_side_configurations",4096)]:
        close(value,summary[field],field)
    return {"status":"PASS","common_genes":len(ga),"cohort_configurations":[32,128],
            "independently_reconstructed_joint_null_values":4096,
            "observed_spearman_rho":observed_rho,"observed_sign_concordance":observed_agree,
            "exact_two_sided_rho_p":rho_p,"exact_two_sided_sign_concordance_p":sign_p,
            "scope":"Exported original-estimator refit coefficients; does not certify raw alignments or anatomical side"}


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report",type=Path)
    args=parser.parse_args()
    result=verify()
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
