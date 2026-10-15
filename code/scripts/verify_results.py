#!/usr/bin/env python3
"""Independently verify retained results with Python's standard library.

This reads tables only: it neither downloads data nor executes the R generators.
It checks schema/probability bounds, BH families, signed scores, exact tests,
diagnosis permutations, sensitivity projections, discovery selection and the
complete patient-side concordance null. This verifier operates on retained
tables; raw-source reconstruction and RNA alignment have separate evidence.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import math
import re
import statistics as st
from collections import defaultdict
from pathlib import Path
from project_paths import resolve_analysis_input

ROOT = Path(__file__).resolve().parents[1]
CHECKS = []
MISSING = {"", "NA", "NaN", "nan", "not retained/detectable"}


def read_table(path):
    # Internal directory scans and synthetic schema tests may supply an already
    # resolved path.  Only repository-relative names participate in the alias map.
    path = Path(path) if Path(path).is_absolute() else resolve_analysis_input(path, ROOT)
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise AssertionError(f"Missing/duplicate column names: {path}")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise AssertionError(f"Malformed TSV row: {path}")
    return rows


def close(actual, expected, context, tolerance=2e-10):
    actual, expected = float(actual), float(expected)
    if not (math.isfinite(actual) and math.isfinite(expected)):
        raise AssertionError(f"Nonfinite numeric comparison: {context}")
    if abs(actual - expected) > tolerance * max(1.0, abs(expected)):
        raise AssertionError(f"{context}: stored={actual:.17g}, independent={expected:.17g}")


def bh(p):
    n = len(p)
    order = sorted(range(n), key=p.__getitem__)
    result = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, p[index] * n / rank)
        result[index] = running
    return result


def _beta_fraction(a, b, x):
    """Continued fraction for the regularized incomplete beta function."""
    small = 1e-300
    qab, qap, qam = a + b, a + 1, a - 1
    c = 1.0
    d = 1 - qab * x / qap
    d = 1 / (d if abs(d) >= small else small)
    h = d
    for m in range(1, 401):
        aa = m * (b - m) * x / ((qam + 2 * m) * (a + 2 * m))
        d = 1 + aa * d
        c = 1 + aa / c
        d = 1 / (d if abs(d) >= small else small)
        c = c if abs(c) >= small else small
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + 2 * m) * (qap + 2 * m))
        d = 1 + aa * d
        c = 1 + aa / c
        d = 1 / (d if abs(d) >= small else small)
        c = c if abs(c) >= small else small
        delta = d * c
        h *= delta
        if abs(delta - 1) < 3e-14:
            return h
    raise ArithmeticError("Incomplete-beta continued fraction did not converge")


def beta_regularized(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    factor = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                      + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1) / (a + b + 2):
        return factor * _beta_fraction(a, b, x) / a
    return 1 - factor * _beta_fraction(b, a, 1 - x) / b


def t_two_sided(value, df):
    return beta_regularized(df / 2, 0.5, df / (df + value * value))


def t_critical(df):
    low, high = 0.0, 64.0
    for _ in range(80):
        middle = (low + high) / 2
        if t_two_sided(middle, df) > 0.05:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def exact_side(d):
    observed = abs(st.mean(d))
    null = [math.fsum(x * s for x, s in zip(d, signs)) / len(d)
            for signs in itertools.product((-1, 1), repeat=len(d))]
    return sum(abs(x) >= observed - 1e-12 for x in null) / len(null)


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + 1 + j) / 2
        for k in order[i:j]:
            result[k] = rank
        i = j
    return result


def rank_unit(values):
    ranked = ranks(values)
    mean = st.mean(ranked)
    centered = [x - mean for x in ranked]
    norm = math.sqrt(math.fsum(x * x for x in centered))
    if norm == 0:
        raise AssertionError("Undefined correlation for constant gene-effect vector")
    return [x / norm for x in centered]


def exact_wilcoxon(d):
    if not all(d) or len(set(abs(x) for x in d)) != len(d):
        raise AssertionError("This verifier expects no zero/tied paired ranks")
    rank = ranks([abs(x) for x in d])
    observed = sum(r for r, x in zip(rank, d) if x > 0)
    total = sum(rank)
    threshold = abs(observed - total / 2)
    sums = [sum(r for r, use in zip(rank, choice) if use)
            for choice in itertools.product((False, True), repeat=len(d))]
    p = sum(abs(v - total / 2) >= threshold - 1e-12 for v in sums) / len(sums)
    return observed, p


def check_summary(d, row, context):
    n = len(d)
    mean, sd = st.mean(d), st.stdev(d)
    se = sd / math.sqrt(n)
    names = {"n_pairs": n, "n_patients": n, "mean_paired_difference": mean,
             "mean_biopsy_minus_culture": mean, "median_paired_difference": st.median(d),
             "paired_difference_SD": sd, "sd_paired_difference": sd,
             "SE": se, "se": se, "cohen_dz": mean / sd,
             "CI95_lower": mean - t_critical(n - 1) * se,
             "CI95_upper": mean + t_critical(n - 1) * se,
             "positive_pairs": sum(x > 0 for x in d),
             "exact_sign_flip_p": exact_side(d),
             "exact_two_sided_sign_flip_p": exact_side(d),
             "t_test_p": t_two_sided(mean / se, n - 1),
             "paired_t_p": t_two_sided(mean / se, n - 1)}
    for column, value in names.items():
        if column in row:
            close(row[column], value, f"{context}/{column}")
    if "exact_wilcoxon_p" in row:
        value, p = exact_wilcoxon(d)
        close(row["wilcoxon_V"], value, context + "/wilcoxon_V")
        close(row["exact_wilcoxon_p"], p, context + "/exact_wilcoxon_p")


def check_pair_map(info, sample_column):
    samples = [r[sample_column] for r in info]
    if len(set(samples)) != len(samples):
        raise AssertionError("Duplicate sample in pairing metadata")
    grouped = defaultdict(dict)
    for row in info:
        side = row["side"]
        if side not in {"concave", "convex"} or side in grouped[row["patient_id"]]:
            raise AssertionError("Missing/duplicate patient-side metadata")
        grouped[row["patient_id"]][side] = row[sample_column]
    if any(set(r) != {"concave", "convex"} for r in grouped.values()):
        raise AssertionError("Incomplete patient pair")
    return dict(sorted(grouped.items()))


def matrix(path):
    rows = read_table(path)
    columns = list(rows[0])[1:]
    result = {}
    for row in rows:
        gene = row[next(iter(row))]
        if gene in result:
            raise AssertionError(f"Duplicate matrix gene in {path}: {gene}")
        values = [float(row[c]) for c in columns]
        if not all(math.isfinite(x) for x in values):
            raise AssertionError(f"Nonfinite expression in {path}: {gene}")
        result[gene] = values
    return columns, result


def standardized(values):
    mean = st.mean(values)
    sd = st.stdev(values)
    return [(x - mean) / sd for x in values] if sd else [0.0] * len(values)


def scores(expression, weights):
    used = {g: w for g, w in weights.items() if g in expression and w != 0}
    denominator = math.fsum(abs(w) for w in used.values())
    n = len(next(iter(expression.values())))
    values = [0.0] * n
    for gene, weight in used.items():
        z = standardized(expression[gene])
        for i in range(n):
            values[i] += weight * z[i] / denominator
    return values, len(used)


def paired(values, columns, info, sample_column):
    values = dict(zip(columns, values))
    pairs = check_pair_map(info, sample_column)
    if set(values) != {v for pair in pairs.values() for v in pair.values()}:
        raise AssertionError("Expression libraries differ from pairing metadata")
    return {p: values[r["concave"]] - values[r["convex"]] for p, r in pairs.items()}


def check_bh(rows, pcolumn, qcolumn, group=()):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[c] for c in group)].append(row)
    for key, family in groups.items():
        expected = bh([float(r[pcolumn]) for r in family])
        for row, value in zip(family, expected):
            close(row[qcolumn], value, f"BH/{key}/{pcolumn}/{qcolumn}")


def schema_and_fdr():
    count, nrows = 0, 0
    directories = [ROOT / "results/tables", ROOT / "results/figures/source_data",
                   ROOT / "analysis/tables", ROOT / "metadata/audit", ROOT / "analysis/bulk_gate_inputs",
                   ROOT / "results/tables/publication"]
    paths = list(ROOT.glob("analysis/*.tsv"))
    paths += [path for directory in directories for path in directory.rglob("*")]
    for path in sorted(set(paths)):
        if not (path.name.endswith(".tsv") or path.name.endswith(".tsv.gz")):
            continue
        rows = read_table(path.relative_to(ROOT))
        count += 1
        nrows += len(rows)
        for row in rows:
            for column, value in row.items():
                if re.search(r"(^PValue|^P\.Value$|^pvalue$|^p_value$|^FDR|_FDR$|^adj\.P\.Val$|_p$|_p_value$)", column) and value not in MISSING:
                    number = float(value)
                    if not math.isfinite(number) or not 0 <= number <= 1:
                        raise AssertionError(f"Invalid probability: {path}/{column}")
        original_full_gene_table = path in {
            ROOT / "analysis/bulk_gate_inputs/original_A_gene_results.tsv.gz",
            ROOT / "analysis/bulk_gate_inputs/original_B_gene_results.tsv.gz",
        }
        if rows and (path.is_relative_to(ROOT / "results/tables") or original_full_gene_table) and {"PValue", "FDR"} <= set(rows[0]):
            group = tuple(c for c in ("cell_type", "contrast") if c in rows[0])
            check_bh(rows, "PValue", "FDR", group)
        if rows and path.is_relative_to(ROOT / "results/tables") and {"p_value", "FDR"} <= set(rows[0]):
            check_bh(rows, "p_value", "FDR")
        if rows and path.is_relative_to(ROOT / "results/tables") and {"t_test_p", "FDR_t_test"} <= set(rows[0]):
            check_bh(rows, "t_test_p", "FDR_t_test")
        if rows and path.is_relative_to(ROOT / "results/tables") and {"exact_sign_flip_p", "FDR_exact"} <= set(rows[0]):
            check_bh(rows, "exact_sign_flip_p", "FDR_exact")
        if rows and path.is_relative_to(ROOT / "results/tables") and {"exact_sign_flip_p", "FDR_exact_within_cohort"} <= set(rows[0]):
            check_bh(rows, "exact_sign_flip_p", "FDR_exact_within_cohort", ("cohort",))
        if rows and path.is_relative_to(ROOT / "results/tables") and {"exact_patient_label_permutation_p", "FDR_interaction"} <= set(rows[0]):
            check_bh(rows, "exact_patient_label_permutation_p", "FDR_interaction")
        if rows and path.is_relative_to(ROOT / "results/tables") and {"meta_p", "meta_FDR"} <= set(rows[0]):
            check_bh(rows, "meta_p", "meta_FDR")
    CHECKS.append({"check": "TSV schemas, probability bounds and original BH families", "files": count, "rows": nrows})



def other_retained_results():
    tables = sorted((ROOT / "results/tables").rglob("*.tsv"))
    paired_count = summary_count = 0
    for path in tables:
        rows = read_table(path)
        if not rows or "paired_difference" not in rows[0]:
            continue
        columns = set(rows[0])
        if {"score.concave", "score.convex"} <= columns:
            first, second = "score.concave", "score.convex"
        elif {"score.biopsy", "score.culture"} <= columns:
            first, second = "score.biopsy", "score.culture"
        else:
            continue
        group = tuple(c for c in ("cell_type", "gene_set", "cohort") if c in columns)
        groups = defaultdict(list)
        for row in rows:
            close(row["paired_difference"], float(row[first]) - float(row[second]), str(path) + "/paired difference")
            groups[tuple(row[c] for c in group)].append(row)
        for family in groups.values():
            assert len({r["patient_id"] for r in family}) == len(family), f"Duplicate patient: {path}"
        paired_count += len(rows)
        candidates = [path.with_name(path.name.replace("paired_effects.tsv", "paired_tests.tsv")),
                      path.with_name(path.name.replace("paired_effects.tsv", "test.tsv")),
                      path.with_name(path.name.replace("joint_paired_effects.tsv", "within_cohort_tests.tsv"))]
        for test_path in dict.fromkeys(candidates):
            if test_path == path or not test_path.is_file():
                continue
            for row in read_table(test_path):
                key = tuple(row[c] for c in group)
                check_summary([float(r["paired_difference"]) for r in groups[key]], row, str(test_path))
                summary_count += 1
    # Descriptive two-donor pathway summaries, including all signs and means.
    for row in read_table("results/tables/programs/scrna_pathway_effect_summary.tsv"):
        a, b = (float(row[c]) for c in ("paired_difference.H01045", "paired_difference.H66206"))
        close(row["mean"], (a + b) / 2, "two-donor pathway mean")
        assert (row["concordant"] in {"TRUE", "1"}) == (a * b > 0)
    composition = read_table("results/tables/composition/cell_type_counts_and_proportions.tsv")
    qc = {row["sample_id"]: row for row in read_table("results/tables/qc/scrna_qc_summary.tsv")}
    samples = defaultdict(list)
    for row in composition:
        samples[row["sample_id"]].append(row)
        close(row["proportion"], int(row["Freq"]) / int(row["sample_total"]), "cell fraction")
    assert sum(int(row["final_singlets"]) for row in qc.values()) == 14846
    for sample, rows in samples.items():
        assert sum(int(r["Freq"]) for r in rows) == int(qc[sample]["final_singlets"])
        assert len({r["cell_type"] for r in rows}) == 9
        assert qc[sample]["ambient_status"] == "success"
        assert int(qc[sample]["pass_qc_metrics"]) - int(qc[sample]["doublets"]) == int(qc[sample]["final_singlets"])
    composition_lookup = {(r["sample_id"], r["cell_type"]): r for r in composition}
    for row in read_table("results/tables/composition/paired_composition_effects.tsv"):
        a = float(composition_lookup[(row["patient_id"] + "__concave", row["cell_type"])]["proportion"])
        b = float(composition_lookup[(row["patient_id"] + "__convex", row["cell_type"])]["proportion"])
        close(row["concave_minus_convex"], a - b, "complete-grid composition difference")
    for path in (ROOT / "results/tables").rglob("deconvolution_two_methods.tsv"):
        groups = defaultdict(list)
        for row in read_table(path):
            value = float(row["proportion"])
            assert math.isfinite(value) and 0 <= value <= 1
            groups[(row["method"], row["sample_id"])].append(value)
        for values in groups.values():
            close(sum(values), 1, "normalized descriptive NNLS weights")
    for row in read_table("results/tables/robustness/pseudobulk_vs_cell_mixed_model.tsv"):
        a, b = float(row["pseudobulk_effect"]), float(row["cell_mixed_effect"])
        assert (row["sign_concordant"] == "TRUE") == ((a > 0) - (a < 0) == (b > 0) - (b < 0))
    embedding = read_table("results/tables/robustness/integration_embedding_sensitivity.tsv")
    for row in read_table("results/tables/robustness/integration_embedding_sensitivity_summary.tsv"):
        for column in ["side_semipartial_r2", "sample_r2"]:
            close(row[column], st.mean(float(r[column]) for r in embedding if r["method"] == row["method"]), "embedding summary")
    for row in read_table("results/tables/communication/cross_method_donor_aware_edges.tsv"):
        retained = all(row[c] in {"1", "TRUE"} for c in ["donor_concordant.natmi", "donor_concordant.cellphonedb",
                        "method_direction_concordant", "pseudobulk_matches_methods"])
        assert (row["retained"] == "TRUE") == retained
    # Independently derive signed Stouffer statistics on the actual common-gene set.
    original_a = {r["gene"]: r for r in read_table("results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz")}
    original_b = {r["gene"]: r for r in read_table("results/tables/bulk_gse311180_requant/GSE311180_paired_edgeR_full.tsv.gz")}
    nd = st.NormalDist()
    for row in read_table("results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz"):
        a, b = original_a[row["gene"]], original_b[row["gene"]]
        for source, suffix in [(a, "GSE254300"), (b, "GSE311180_requant")]:
            for column in ["logFC", "PValue", "FDR"]:
                close(row[column + "." + suffix], source[column], "meta input gene coefficient/" + suffix)
        za = math.copysign(-nd.inv_cdf(max(float(a["PValue"]) / 2, 1e-300)), float(a["logFC"])) if float(a["logFC"]) else 0
        zb = math.copysign(-nd.inv_cdf(max(float(b["PValue"]) / 2, 1e-300)), float(b["logFC"])) if float(b["logFC"]) else 0
        z = (math.sqrt(5) * za + math.sqrt(7) * zb) / math.sqrt(12)
        close(row["meta_z"], z, "signed Stouffer Z")
        close(row["meta_p"], math.erfc(abs(z) / math.sqrt(2)), "signed Stouffer P")
    CHECKS.append({"check": "All stored score-pair differences, available pathway/axis tests, QC counts, descriptive fractions, model signs, communication filters and signed Stouffer arithmetic",
                   "paired_rows": paired_count, "paired_test_rows": summary_count})


def verify(full_null=True):
    CHECKS.clear()
    schema_and_fdr()
    other_retained_results()
    signature = read_table("results/tables/bulk/frozen_multicellular_signature.tsv")
    grouped = defaultdict(list)
    for row in signature:
        grouped[row["gene"]].append(float(row["weight"]))
    weights = {g: st.mean(w) for g, w in grouped.items() if st.mean(w) != 0}
    assert len(signature) == 150 and len(weights) == 144 and len(grouped) == 146
    assert len({(r["gene"], r["cell_type"]) for r in signature}) == 150
    assert len({r["cell_type"] for r in signature}) == 5
    ca, ma = matrix("results/tables/bulk/GSE254300_AIS_label_blind_projection_logCPM.tsv.gz")
    cb, mb = matrix("results/tables/bulk_gse311180_requant/GSE311180_label_blind_projection_logCPM.tsv.gz")
    ia = [r for r in read_table("metadata/audit/bulk_samples.tsv") if r["cohort"] == "AIS"]
    ib = [r for r in read_table("metadata/audit/gse311180_samples.tsv") if r["sample_id"] in cb]
    sa, na = scores(ma, weights)
    sb, nb = scores(mb, weights)
    assert (na, nb) == (122, 126)
    ea = paired(sa, ca, ia, "processed_sample_id")
    eb = paired(sb, cb, ib, "sample_id")
    for path, values in [("results/tables/bulk/multicellular_axis_scores.tsv", dict(zip(ca, sa))),
                         ("results/tables/bulk_gse311180_requant/frozen_axis_scores.tsv", dict(zip(cb, sb)))]:
        for row in read_table(path):
            close(row["score"], values[row["sample_id"]], "original score/" + row["sample_id"])
    ordered_means = sorted(st.mean(x) for x in ma.values())
    position = .99 * (len(ordered_means) - 1)
    index = int(position)
    cutoff = ordered_means[index] + (position - index) * (ordered_means[min(index + 1, len(ordered_means) - 1)] - ordered_means[index])
    trimmed = {g: w for g, w in weights.items() if g in ma and st.mean(ma[g]) < cutoff}
    high_scores, high_count = scores(ma, trimmed)
    high_mean = st.mean(paired(high_scores, ca, ia, "processed_sample_id").values())
    for row in read_table("results/tables/robustness/bulk_axis_high_expression_sensitivity.tsv"):
        count, mean = (na, st.mean(ea.values())) if row["analysis"] == "all_signature_genes" else (high_count, high_mean)
        close(row["n_genes"], count, "high-expression sensitivity actual nonzero gene count")
        close(row["mean_paired_difference"], mean, "high-expression sensitivity projection")
    ordered_b_means = sorted(st.mean(x) for x in mb.values())
    position_b = .99 * (len(ordered_b_means) - 1)
    index_b = int(position_b)
    cutoff_b = ordered_b_means[index_b] + (position_b - index_b) * (
        ordered_b_means[min(index_b + 1, len(ordered_b_means) - 1)] - ordered_b_means[index_b])
    trimmed_b = {g: w for g, w in weights.items() if g in mb and st.mean(mb[g]) < cutoff_b}
    high_scores_b, high_count_b = scores(mb, trimmed_b)
    high_effects_b = list(paired(high_scores_b, cb, ib, "sample_id").values())
    high_mean_b = st.mean(high_effects_b)
    high_se_b = st.stdev(high_effects_b) / math.sqrt(len(high_effects_b))
    high_rows_b = read_table("results/tables/bulk_gse311180_requant/frozen_axis_remove_top1pct.tsv")
    if len(high_rows_b) != 1:
        raise AssertionError("Ambiguous seven-patient high-expression sensitivity")
    for field, value in [("n_signature_genes", nb), ("n_retained_genes", high_count_b),
                         ("mean_paired_difference", high_mean_b),
                         ("positive_pairs", sum(x > 0 for x in high_effects_b)),
                         ("t_test_p", t_two_sided(high_mean_b / high_se_b, len(high_effects_b) - 1))]:
        close(high_rows_b[0][field], value, "seven-patient high-expression sensitivity/" + field)
    pooled = read_table("results/tables/meta_validation_requant/frozen_axis_12_patient_effects.tsv")
    for row in pooled:
        close(row["paired_difference"], (ea if row["cohort"] == "GSE254300" else eb)[row["patient_id"]], "pooled paired score")
    d = [float(r["paired_difference"]) for r in pooled]
    for row in read_table("analysis/bulk_patient_effect_statistics.tsv"):
        values = list(ea.values()) if row["analysis"] == "GSE254300" else list(eb.values()) if row["analysis"] == "GSE311180_requant" else d
        check_summary(values, row, row["analysis"])
    for row in read_table("analysis/bulk_pooled_leave_one_patient_out.tsv"):
        values = [float(x["paired_difference"]) for x in pooled if (x["cohort"], x["patient_id"]) != (row["omitted_cohort"], row["omitted_patient"])]
        check_summary(values, row, "pooled leave-one-out")
    CHECKS.append({"check": "144 weights, 122/126 projections, pooled summaries and all 12 omissions", "status": "pass"})

    omitted = {r["patient_id"] for r in read_table("analysis/audit_gse311180_rna_expected_pairs.tsv")
               if r["interpretation"] != "consistent_RNA_profiles"}
    assert omitted == {"STB-084", "STB-100"}
    kept = [i for i, sample in enumerate(cb) if next(r["patient_id"] for r in ib if r["sample_id"] == sample) not in omitted]
    columns_kept = [cb[i] for i in kept]
    matrix_kept = {gene: [value[i] for i in kept] for gene, value in mb.items()}
    info_kept = [r for r in ib if r["patient_id"] not in omitted]
    rez_scores, _ = scores(matrix_kept, weights)
    rez_effects = paired(rez_scores, columns_kept, info_kept, "sample_id")
    locked_effects = {p: value for p, value in eb.items() if p not in omitted}
    for row in read_table("analysis/bulk_RNA_pairing_gate_exclusion.tsv"):
        selected = locked_effects if "locked" in row["analysis"] else rez_effects
        values = list(selected.values())
        if "pooled" in row["analysis"]:
            values = list(ea.values()) + values
        check_summary(values, row, "RNA gate exclusion/" + row["analysis"])
        assert set(row["omitted_patients"].split(";")) == omitted
    for row in read_table("analysis/bulk_RNA_pairing_gate_exclusion_patient_effects.tsv"):
        selected = locked_effects if row["scaling"] == "locked" else rez_effects
        close(row["paired_difference"], selected[row["patient_id"]], "RNA excluded-pair score")
    CHECKS.append({"check": "Fixed-gate RNA pair exclusions, locked/restandardized projections and exact tests", "status": "pass"})

    # Audit subset weights against the actual gene-z matrices, rather than only
    # reproducing a P value from an already generated paired-effect file.
    projected = read_table("analysis/signature_bulk_sensitivity_projected_weights.tsv")
    subsets = defaultdict(dict)
    for row in projected:
        key = (row["analysis"], row["cohort"])
        if row["gene"] in subsets[key]:
            raise AssertionError("Duplicated sensitivity gene")
        subsets[key][row["gene"]] = float(row["weight"])
    sensitivity = {}
    for key, w in subsets.items():
        c, m, info, samplecolumn = (ca, ma, ia, "processed_sample_id") if key[1] == "GSE254300_AIS" else (cb, mb, ib, "sample_id")
        values, n = scores(m, w)
        sensitivity[key] = paired(values, c, info, samplecolumn)
    for row in read_table("analysis/signature_bulk_sensitivity_tests.tsv"):
        key = (row["analysis"], row["cohort"])
        if row["cohort"] == "pooled_12_AIS":
            values = list(sensitivity[(row["analysis"], "GSE254300_AIS")].values()) + list(sensitivity[(row["analysis"], "GSE311180_AIS")].values())
        else:
            values = list(sensitivity[key].values())
        check_summary(values, row, "sensitivity/" + str(key))
    subrows = [r for r in read_table("analysis/signature_bulk_sensitivity_tests.tsv") if r["analysis"].startswith("lineage_only__")]
    check_bh(subrows, "exact_two_sided_sign_flip_p", "FDR_five_lineage_subscores", ("cohort",))
    CHECKS.append({"check": "All marker/lineage sensitivity projections and five-lineage BH families", "variants": len(subsets)})

    # Reselection uses each whole donor vector, never independent gene shuffles.
    de = read_table("results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz")
    cols = ["paired_logCPM_difference.H01045", "paired_logCPM_difference.H66206"]
    supplied = read_table("analysis/signature_pipeline_flip_reselected_signatures.tsv")
    for signs in itertools.product((1, -1), repeat=2):
        label = ",".join(map(str, signs))
        by_lineage = defaultdict(list)
        for row in de:
            a, b = (float(row[c]) * s for c, s in zip(cols, signs))
            if a * b > 0 and not re.match(r"^(MT-|RPS|RPL|HBA|HBB)", row["gene"]):
                by_lineage[row["cell_type"]].append((row["gene"], (a + b) / 2))
        selected = {(g, ct): (i + 1, math.copysign(1, effect), effect)
                    for ct, genes in by_lineage.items()
                    for i, (g, effect) in enumerate(sorted(genes, key=lambda x: -abs(x[1]))[:30])}
        rows = [r for r in supplied if r["configuration"] == label]
        assert set(selected) == {(r["gene"], r["cell_type"]) for r in rows}
        for row in rows:
            rank, weight, effect = selected[(row["gene"], row["cell_type"])]
            close(row["rank_within_lineage"], rank, "discovery rank")
            close(row["weight"], weight, "discovery sign")
            close(row["mean_donor_effect"], effect, "discovery effect")
    CHECKS.append({"check": "All four whole-donor discovery reselections", "status": "pass"})

    joint = read_table("results/tables/disease_control_gse254300_cs/frozen_axis_joint_paired_effects.tsv")
    a = [float(r["paired_difference"]) for r in joint if r["cohort"] == "AIS"]
    b = [float(r["paired_difference"]) for r in joint if r["cohort"] == "CS"]
    both = a + b
    observed = st.mean(a) - st.mean(b)
    null = []
    for chosen in itertools.combinations(range(10), 5):
        chosen = set(chosen)
        null.append(st.mean([x for i, x in enumerate(both) if i in chosen]) - st.mean([x for i, x in enumerate(both) if i not in chosen]))
    row = read_table("analysis/bulk_interaction_uncertainty.tsv")[0]
    close(row["difference"], observed, "diagnosis interaction")
    close(row["exact_patient_label_permutation_p"], sum(abs(x) >= abs(observed) - 1e-12 for x in null) / 252, "diagnosis exact P")
    va, vb = st.variance(a) / 5, st.variance(b) / 5
    se = math.sqrt(va + vb)
    df = (va + vb)**2 / (va**2 / 4 + vb**2 / 4)
    close(row["welch_df"], df, "Welch df")
    close(row["CI95_lower"], observed - t_critical(df) * se, "Welch lower")
    close(row["CI95_upper"], observed + t_critical(df) * se, "Welch upper")
    CHECKS.append({"check": "252 diagnosis permutations and independently computed Welch interval", "status": "pass"})

    # Empirical random-signature tails are checked directly from all 10,000 draws.
    random = read_table("analysis/bulk_pooled_matched_random_signature_null.tsv")
    for row in read_table("analysis/bulk_pooled_matched_random_signature_summary.tsv"):
        values = [float(x[row["analysis"]]) for x in random]
        observed = float(row["observed_mean"])
        close(row["empirical_two_sided_p"], (1 + sum(abs(x) >= abs(observed) for x in values)) / (1 + len(values)), "joint random tail")
        close(row["null_mean"], st.mean(values), "joint null mean")
        close(row["null_SD"], st.stdev(values), "joint null SD")
    CHECKS.append({"check": "10,000-draw pooled/cohort empirical random-signature tails", "status": "pass"})

    rep = read_table("results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz")
    common = [r["gene"] for r in rep]
    xa, ta = matrix("results/tables/bulk/GSE254300_AIS_TMM_logCPM.tsv.gz")
    xb, tb = matrix("results/tables/bulk_gse311180_requant/GSE311180_TMM_logCPM.tsv.gz")
    pa, pb = check_pair_map(ia, "processed_sample_id"), check_pair_map(ib, "sample_id")
    def differences(columns, expressions, pairs):
        positions = {s: i for i, s in enumerate(columns)}
        return [[expressions[g][positions[r["concave"]]] - expressions[g][positions[r["convex"]]] for r in pairs.values()] for g in common]
    da, db = differences(xa, ta, pa), differences(xb, tb, pb)
    means_a, means_b = [st.mean(x) for x in da], [st.mean(x) for x in db]
    ra, rb = rank_unit(means_a), rank_unit(means_b)
    rho = math.fsum(x * y for x, y in zip(ra, rb))
    sign = lambda x: (x > 0) - (x < 0)
    agreement = sum(sign(x) == sign(y) for x, y in zip(means_a, means_b)) / len(common)
    row = read_table("analysis/bulk_gene_concordance_patient_permutation_summary.tsv")[0]
    close(row["observed_spearman_rho"], rho, "concordance observed rho")
    close(row["observed_sign_concordance"], agreement, "concordance observed signs")
    if full_null:
        def perm_vectors(d, n):
            # R expand.grid reverses the order of itertools.product columns.
            vectors = [[math.fsum(x * s for x, s in zip(values, reversed(signs))) / n for values in d]
                       for signs in itertools.product((-1, 1), repeat=n)]
            return [rank_unit(x) for x in vectors], [[sign(v) for v in x] for x in vectors]
        ranked_a, signs_a = perm_vectors(da, 5)
        ranked_b, signs_b = perm_vectors(db, 7)
        null_rows = read_table("analysis/bulk_gene_concordance_patient_permutation_null.tsv")
        rho_tail, sign_tail = 0, 0
        for nr in null_rows:
            i = int(nr["GSE254300_configuration"]) - 1
            j = int(nr["GSE311180_configuration"]) - 1
            value = math.fsum(x * y for x, y in zip(ranked_a[i], ranked_b[j]))
            agree = sum(x == y for x, y in zip(signs_a[i], signs_b[j])) / len(common)
            close(nr["spearman_rho"], value, "complete rho null")
            close(nr["sign_concordance"], agree, "complete sign null")
            rho_tail += abs(value) >= abs(rho) - 1e-12
            sign_tail += abs(agree - .5) >= abs(agreement - .5) - 1e-12
        close(row["exact_two_sided_rho_p"], rho_tail / 4096, "rho exact P")
        close(row["exact_two_sided_sign_concordance_p"], sign_tail / 4096, "sign exact P")
    CHECKS.append({"check": "Independent TMM correlation/sign computation and complete null" if full_null else "Independent observed TMM correlation/sign computation", "configurations": 4096 if full_null else 1})
    return {"status": "PASS", "checks": CHECKS,
            "limits": ["Raw-data DE/QC/cell-expression pipelines cannot be verified from retained tables alone",
                       "This verification does not establish specimen identity or biological validity",
                       "No prospective freeze or public repository identifier is certified"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--skip-concordance-null", action="store_true")
    args = parser.parse_args()
    report = verify(not args.skip_concordance_null)
    text = json.dumps(report, indent=2) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(text, encoding="utf-8")
    print(text)
