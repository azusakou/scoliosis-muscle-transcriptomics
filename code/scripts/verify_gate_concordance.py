#!/usr/bin/env python3
"""Reconstruct the fixed RNA-gate concordance nulls and candidate summaries.

All 1,024 patient-side configurations are checked for each normalization.
The R generator separately refits the original and retained count models.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import statistics
from pathlib import Path
from verify_edgeR_concordance_null import read_coefficients, rank_unit, close, table
from verify_results import bh, t_critical, exact_side

ROOT = Path(__file__).resolve().parents[1]
NORMALIZATIONS = ("locked_original_TMM", "retained_library_reTMM")
REQUIRED_GATE_INPUTS = frozenset({
    "original_A_gene_results.tsv.gz", "original_B_gene_results.tsv.gz",
    "original_B_QC.tsv", "A32_logFC.tsv.gz", "A32_grid.tsv",
    "RNA_expected_pairs.tsv", "GSE254300_tested_counts.tsv.gz",
    "GSE311180_tested_counts.tsv.gz", "samples.tsv", "common_gene_universe.tsv",
    "published_retained_TMM_factors.tsv",
})


def verify_input_bindings(inputs):
    """Every used input, including fixed factors, must be explicitly SHA-bound."""
    bindings = list(table(inputs / "source_bindings.tsv")); seen = set()
    for row in bindings:
        name = row["artifact"]
        if not name or name in seen or Path(name).name != name:
            raise ValueError("Invalid gate input path")
        seen.add(name)
        if not re.fullmatch(r"[0-9a-f]{64}", row["SHA256"]):
            raise ValueError("Invalid gate input SHA256: " + name)
        if hashlib.sha256((inputs/name).read_bytes()).hexdigest() != row["SHA256"]:
            raise ValueError("Gate input hash mismatch: " + name)
    if not REQUIRED_GATE_INPUTS <= seen:
        raise ValueError("Unbound required gate inputs: " + ", ".join(sorted(REQUIRED_GATE_INPUTS-seen)))
    return bindings


def unique_samples(records, label, required):
    indexed = {}
    for row in records:
        if not required <= set(row):
            raise ValueError("Missing normalization source fields: " + label)
        sample = row["sample_id"]
        if not sample or sample in indexed or not row["patient_id"] or row["side"] not in ("concave", "convex"):
            raise ValueError("Invalid/duplicate normalization sample ID: " + label)
        indexed[sample] = row
    if not indexed:
        raise ValueError("Empty normalization source: " + label)
    return indexed


def positive_finite(value, label):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("Nonpositive/nonfinite normalization value: " + label) from exc
    if not math.isfinite(number) or number <= 0:
        raise ValueError("Nonpositive/nonfinite normalization value: " + label)
    return number


def verify_normalization_join(qc, published, reported, metadata, original_map):
    """Verify fixed model parameters by identity, not row order or fresh TMM."""
    identity = {"sample_id", "patient_id", "side"}
    factor_fields = {"locked_original_norm_factor", "retained_library_TMM_factor", "retained_tested_library_size"}
    qc_index = unique_samples(qc, "original QC", identity | {"normalization_factor"})
    published_index = unique_samples(published, "published retained factors", identity | factor_fields)
    reported_index = unique_samples(reported, "reported factors", identity | factor_fields)
    metadata_index = unique_samples(metadata, "subset metadata", identity)
    original_index = unique_samples(original_map, "original sample map", identity)
    if set(qc_index) != set(original_index):
        raise ValueError("Original QC normalization sample coverage changed")
    expected = set(metadata_index)
    if set(published_index) != expected or set(reported_index) != expected or not expected <= set(original_index):
        raise ValueError("Retained normalization sample coverage changed")
    for sample, row in qc_index.items():
        if any(row[field] != original_index[sample][field] for field in identity):
            raise ValueError("Original QC normalization identity mismatch: " + sample)
        positive_finite(row["normalization_factor"], "original QC " + sample)
    for sample, info in metadata_index.items():
        if any(info[field] != original_index[sample][field] for field in identity):
            raise ValueError("Subset normalization identity mismatch: " + sample)
        p, r = published_index[sample], reported_index[sample]
        for label, row in (("published", p), ("reported", r)):
            if any(row[field] != info[field] for field in identity):
                raise ValueError(label + " normalization identity mismatch: " + sample)
            for field in factor_fields:
                positive_finite(row[field], label + " " + sample + " " + field)
        close(r["locked_original_norm_factor"], qc_index[sample]["normalization_factor"], "locked factor QC join " + sample)
        close(p["locked_original_norm_factor"], qc_index[sample]["normalization_factor"], "published locked factor QC join " + sample)
        close(r["retained_library_TMM_factor"], p["retained_library_TMM_factor"], "published retained factor join " + sample)
        close(r["retained_tested_library_size"], p["retained_tested_library_size"], "published tested-library size join " + sample)
    return dict(locked_factors_joined=len(expected), published_retained_factors_joined=len(expected),
                original_QC_sample_ids_verified=len(original_index),
                tolerance="unchanged close arithmetic tolerance 2e-10; native recalculation is diagnostic only")


def exact_grid(records, patient_count):
    names = [k for k in records[0] if k != "configuration"]
    if len(names) != patient_count or len(records) != 2**patient_count:
        raise ValueError("Invalid gate sensitivity grid size")
    for i, row in enumerate(records):
        if int(row["configuration"]) != i+1 or tuple(int(row[k]) for k in names) != tuple(
                1 if (i >> j) & 1 else -1 for j in range(patient_count)):
            raise ValueError("Invalid gate patient-side configuration")


def verify(root=ROOT, require_generated_normalization=False):
    rev = root / "analysis"; inputs = rev / "bulk_gate_inputs"
    bindings = verify_input_bindings(inputs)
    ga, a = read_coefficients(inputs / "A32_logFC.tsv.gz", 32)
    common = [r["gene"] for r in table(inputs / "common_gene_universe.tsv")]
    if ga != common or len(set(common)) != 12995: raise ValueError("Fixed common universe changed")
    exact_grid(list(table(inputs / "A32_grid.tsv")), 5)
    exact_grid(list(table(rev / "bulk_gate_B32_grid.tsv")), 5)
    gate = list(table(inputs / "RNA_expected_pairs.tsv"))
    # The compact release retains one SHA-bound copy of the external decisions;
    # verify_rna_identity.py independently reconstructs them from all-pair data.
    keep = {r["patient_id"] for r in gate if r["interpretation"] == "consistent_RNA_profiles"}
    metadata = list(table(rev / "bulk_gate_subset_metadata.tsv"))
    second_map = [r for r in table(inputs / "samples.tsv") if r["cohort"] == "GSE311180"]
    original_map = [r for r in second_map if r["patient_id"] in keep]
    if metadata != original_map or len(keep) != 5 or len(metadata) != 10: raise ValueError("Gate subset source join")
    published = list(table(inputs / "published_retained_TMM_factors.tsv"))
    reported_path = rev / "bulk_gate_normalization_factors.tsv"
    if require_generated_normalization and not reported_path.is_file():
        raise ValueError("Freshly generated normalization output is required")
    generated = reported_path.is_file()
    if reported_path.is_symlink():
        raise ValueError("Symlink normalization output")
    normalization_join = verify_normalization_join(
        list(table(inputs / "original_B_QC.tsv")),
        published, list(table(reported_path)) if generated else published,
        metadata, second_map)
    normalization_join["verification_mode"] = (
        "fresh generated output joined against separate fixed published input" if generated
        else "fixed published input identity/QC check only; generated output is not shipped in the compact release")
    normalization_join["fresh_generated_output_verified"] = generated
    if not generated:
        normalization_join["published_retained_factors_joined"] = 0
    summaries = list(table(rev / "bulk_gate_genomewide_exact_summary.tsv"))
    indexed = {(r["normalization"], r["endpoint"]): r for r in summaries}
    if len(indexed) != 4 or len(summaries) != 4: raise ValueError("Endpoint family coverage")
    primary_a = [r["logFC"] for r in table(inputs / "original_A_gene_results.tsv.gz")]
    olda = {r["gene"]: r for r in table(inputs / "original_A_gene_results.tsv.gz")}
    for i, gene in enumerate(ga): close(a[-1][i], olda[gene]["logFC"], "original A coefficient")
    ra = [rank_unit(col) for col in a]; signs_a = [[(v > 0)-(v < 0) for v in c] for c in a]
    new_tables = {}; reconstructed = []
    for norm in NORMALIZATIONS:
        gb, b = read_coefficients(rev / ("bulk_gate_B32_logFC_"+norm+".tsv.gz"), 32)
        full = list(table(rev / ("bulk_gate_observed_genes_"+norm+".tsv.gz")))
        original_genes = [r["gene"] for r in table(inputs / "original_B_gene_results.tsv.gz")]
        if gb != ga or [r["gene"] for r in full] != original_genes or len(full) != 17961:
            raise ValueError("Subset refit changed gene universe")
        new_tables[norm] = {r["gene"]: r for r in full}
        for i, gene in enumerate(gb): close(b[-1][i], new_tables[norm][gene]["logFC"], "observed B coefficient")
        full_bh = bh([float(r["PValue"]) for r in full])
        for r, adj in zip(full, full_bh): close(adj, r["FDR"], "full gene-family BH")
        rb = [rank_unit(col) for col in b]; signs_b = [[(v > 0)-(v < 0) for v in c] for c in b]
        null = list(table(rev / ("bulk_gate_exact_null_"+norm+".tsv"))); visited = set(); rhos = []; signs = []
        for row in null:
            i, j = int(row["A_configuration"])-1, int(row["B_configuration"])-1
            if not (0 <= i < 32 and 0 <= j < 32) or (i, j) in visited: raise ValueError("Duplicate/invalid exact null")
            visited.add((i,j))
            rho = math.fsum(x*y for x,y in zip(ra[i], rb[j]))
            agree = sum(x == y for x,y in zip(signs_a[i], signs_b[j]))/len(ga)
            close(rho, row["spearman_rho"], "gate exact rho"); close(agree, row["sign_concordance"], "gate exact signs")
            rhos.append(rho); signs.append(agree)
        if len(visited) != 1024: raise ValueError("Incomplete gate exact null")
        rho = math.fsum(x*y for x,y in zip(ra[-1], rb[-1])); agree = sum(x==y for x,y in zip(signs_a[-1], signs_b[-1]))/len(ga)
        ps = [sum(abs(v) >= abs(rho)-1e-12 for v in rhos)/1024,
              sum(abs(v-.5) >= abs(agree-.5)-1e-12 for v in signs)/1024]
        for endpoint, obs, pv, adj in zip(("spearman_logFC", "sign_concordance"), (rho, agree), ps, bh(ps)):
            row = indexed[norm, endpoint]
            for field, value in (("observed", obs), ("exact_two_sided_p", pv), ("BH_within_two_endpoint_family", adj),
                                 ("n_A_pairs", 5), ("n_B_pairs", 5), ("n_common_genes", 12995), ("n_joint_side_configurations", 1024)):
                close(value, row[field], field)
        reconstructed.append(dict(normalization=norm, rho=rho, sign_concordance=agree, exact_p=ps))
    for r, adj in zip(summaries, bh([float(r["exact_two_sided_p"]) for r in summaries])):
        close(adj, r["BH_full_four_endpoint_sensitivity_family"], "four-endpoint BH")
    candidates = list(table(rev / "bulk_gate_ANGPT1_LEP_fixed_candidate_results.tsv"))
    patient_values = list(table(rev / "bulk_gate_fixed_candidate_patient_effects.tsv"))
    exact_rows = list(table(rev / "bulk_gate_candidate_patient_exact.tsv"))
    exact_lookup = {(r["normalization"], r["gene"]): r for r in exact_rows}
    if len(exact_lookup) != 4 or len(exact_rows) != 4: raise ValueError("Candidate exact-test family coverage")
    if {(r["normalization"],r["gene"]) for r in candidates} != {(n,g) for n in NORMALIZATIONS for g in ("ANGPT1","LEP")} or len(candidates) != 4:
        raise ValueError("Candidate family changed")
    for norm in NORMALIZATIONS:
        selected = [r for r in candidates if r["normalization"] == norm]
        for r, adj in zip(selected, bh([float(r["QL_P"]) for r in selected])):
            original = new_tables[norm][r["gene"]]
            for source, dest in (("logFC", "fitted_shrunk_logFC"), ("F", "QL_F"), ("PValue", "QL_P"), ("FDR", "genome_wide_FDR_original_17961")):
                close(float(original[source]), r[dest], "candidate full-model source")
            close(adj, r["BH_two_fixed_candidates_within_normalization"], "candidate BH2")
            values = [p for p in patient_values if p["normalization"] == norm and p["gene"] == r["gene"]]
            if len(values) != 5 or {p["patient_id"] for p in values} != keep: raise ValueError("Candidate patient identity")
            d = [float(p["paired_logCPM_difference"]) for p in values]; mean = statistics.fmean(d); sd = statistics.stdev(d)
            half = t_critical(4)*sd/math.sqrt(5)
            exact = exact_lookup[norm, r["gene"]]
            for field, value in (("n_pairs",5),("positive_pairs",sum(v>0 for v in d)),
                                 ("paired_logCPM_mean",mean),("exact_two_sided_p",exact_side(d)),
                                 ("minimum_attainable_two_sided_p",2/32)):
                close(value,exact[field],field)
            for field, value in (("patient_paired_logCPM_mean",mean),("patient_paired_logCPM_SD",sd),
                                 ("patient_t_CI95_lower",mean-half),("patient_t_CI95_upper",mean+half),
                                 ("positive_pairs",sum(v>0 for v in d)),("n_pairs",5)):
                close(value,r[field],field)
    for r, adj in zip(candidates,bh([float(r["QL_P"]) for r in candidates])):
        close(adj,r["BH_four_candidate_normalization_sensitivity_family"],"candidate BH4")
    for r, adj in zip(exact_rows,bh([float(r["exact_two_sided_p"]) for r in exact_rows])):
        close(adj,r["BH_four_candidate_normalization_family"],"candidate exact BH4")
    return dict(status="PASS", fixed_genes=12995, B_refit_genes_per_normalization=17961,
                configurations_per_normalization=1024, independently_reconstructed_null_entries=2048,
                normalization_results=reconstructed, source_bound_portable_inputs=len(bindings),
                fixed_normalization_parameter_join=normalization_join,
                fixed_candidates=2, scope="Exported complete refits/metadata/count-input hashes; no signature rescue or anatomical authentication")


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--root",type=Path,default=ROOT); p.add_argument("--report",type=Path)
    p.add_argument("--require-generated-normalization", action="store_true",
                   help="Require a new run output, distinct from the immutable published-factor input")
    args=p.parse_args(); result=verify(args.root, args.require_generated_normalization)
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True); args.report.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
