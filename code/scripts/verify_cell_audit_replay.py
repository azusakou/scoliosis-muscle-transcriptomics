#!/usr/bin/env python3
"""Compare full fresh cell audits and run the two independent cell verifiers.

Counts and labels remain exact. Numerical comparisons use the existing ordinary
arithmetic criteria; this module does not introduce a tolerance or a model fit.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from collections import Counter
from pathlib import Path
import tempfile

import retained_reproduction as reproduction
import verify_biology_cells
import verify_cell_expression


ROOT = Path(__file__).resolve().parents[1]
SCIENCE_OUTPUTS = (
    "signature_requested_gene_availability.tsv",
    "signature_expression_by_library_lineage_before_after.tsv",
    "signature_unexpected_gene_cell_source.tsv.gz",
    "signature_expression_correction_changes.tsv",
    "signature_frozen_gene_expression_all_four_libraries.tsv",
    "signature_positive_cell_QC_trace.tsv.gz",
    "signature_positive_cell_QC_concentration.tsv",
    "signature_biology_target_panel_manifest.tsv",
    "signature_biology_target_gene_expression_all_libraries_lineages.tsv",
    "signature_biology_all_target_cell_source.tsv.gz",
    "signature_biology_bulk_candidate_cell_source.tsv.gz",
    "signature_biology_target_gene_descriptive_within_donor_differences.tsv",
    "signature_biology_retained_cell_expression_lineage_provenance.tsv",
)


# These are deterministic expansions of authentic sparse counts, not extra inputs.
DERIVED_CELLS = {
    "signature_unexpected_gene_cell_source.tsv.gz": verify_cell_expression.GENES,
    "signature_positive_cell_QC_trace.tsv.gz": verify_cell_expression.GENES,
    "signature_biology_all_target_cell_source.tsv.gz": verify_biology_cells.TARGETS,
    "signature_biology_bulk_candidate_cell_source.tsv.gz": {"ANGPT1", "LEP", "ADIPOQ"},
}


def check_sparse_row(row, metadata, counts, available):
    assay, cell, gene = row["assay"], row["cell_id"], row["gene"]
    if (assay, gene) not in available or cell not in metadata:
        raise ValueError("Unavailable gene/assay or unknown cell")
    info = metadata[cell]
    for field in ("sample_id", "patient_id", "side", "cell_type"):
        if row[field] != info[field]:
            raise ValueError("Sparse-source cell label mismatch")
    value = float(row["UMI"])
    if not math.isfinite(value) or not value.is_integer() or value != counts.get((assay, cell, gene), 0):
        raise ValueError("Cell UMI differs from authentic sparse input")
    total = float(info["nCount_" + assay])
    if not math.isfinite(total) or total <= 0 or value > total:
        raise ValueError("Invalid full-assay cell total")
    if "all_gene_total_UMI" in row and float(row["all_gene_total_UMI"]) != total:
        raise ValueError("Changed full-assay cell total")
    normalized = value / total * 10000
    for field, expected in (("cp10k", normalized), ("log1p_cp10k", math.log1p(normalized))):
        if field in row and not math.isclose(float(row[field]), expected,
                rel_tol=reproduction.REL_TOL, abs_tol=reproduction.ABS_TOL):
            raise ValueError("Sparse-source normalization mismatch")
    for field in ("nCount_RNA", "nFeature_RNA", "nCount_SoupX", "nFeature_SoupX"):
        if field in row and float(row[field]) != float(info[field]):
            raise ValueError("Changed exact QC count")
    for field in ("subsets_Mito_percent", "subsets_Ribo_percent", "scDblFinder.score"):
        if field in row and not math.isclose(float(row[field]), float(info[field]),
                rel_tol=reproduction.REL_TOL, abs_tol=reproduction.ABS_TOL):
            raise ValueError("Changed QC source value")
    for field in ("scDblFinder.class", "passes_qc_metrics", "passes_qc", "seurat_clusters", "annotation_confidence"):
        if field in row and row[field] != info[field]:
            raise ValueError("Changed QC source label")
    return assay, cell, gene


def verify_sparse_expansions(root, run):
    rows = verify_biology_cells.rows
    metadata_rows = list(rows(root / "results/tables/annotation/annotated_cell_metadata.tsv.gz"))
    metadata = {r["cell_id"]: r for r in metadata_rows}
    if len(metadata) != len(metadata_rows) or len(metadata) != 14846:
        raise ValueError("Sparse-source metadata coverage")
    available = {(r["assay"], r["gene"]) for r in rows(root / "data/retained/cell_target_availability.tsv")
                 if r["available"] == "TRUE"}
    counts = {}
    for row in rows(root / "data/retained/cell_target_counts.tsv.gz"):
        key = row["assay"], row["cell_id"], row["gene"]
        value = float(row["UMI"])
        if key in counts or key[:1] + key[2:] not in available or key[1] not in metadata or not math.isfinite(value) or value <= 0 or not value.is_integer():
            raise ValueError("Invalid authentic sparse count")
        counts[key] = int(value)
    result = []
    for name, genes in DERIVED_CELLS.items():
        seen = set(); coverage = Counter()
        positive_only = name == "signature_positive_cell_QC_trace.tsv.gz"
        for row in rows(run / "analysis" / name):
            key = check_sparse_row(row, metadata, counts, available)
            if key[2] not in genes or key in seen or (positive_only and counts.get(key, 0) <= 0):
                raise ValueError("Unexpected or duplicate sparse expansion record")
            seen.add(key); coverage[(key[0], key[2])] += 1
        expected = sum(1 for key in counts if key[2] in genes) if positive_only else 2 * len(metadata) * len(genes)
        if len(seen) != expected:
            raise ValueError("Incomplete sparse expansion coverage")
        if not positive_only and (set(coverage) != {(a, g) for a in ("RNA", "SoupX") for g in genes} or
                any(n != len(metadata) for n in coverage.values())):
            raise ValueError("Incomplete cell/assay/gene coverage")
        result.append(dict(file=name, records=len(seen), direct_sparse_input_join="PASS"))
    return result


def assay_manifest_comparison(root: Path, run: Path) -> dict:
    """Compare all scientific fields; require an honest newly declared source."""
    with tempfile.TemporaryDirectory() as directory:
        temporary = Path(directory)
        for origin, name in ((root, "reference"), (run, "fresh")):
            with (origin / "analysis/signature_expression_assay_manifest.tsv").open(newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields = reader.fieldnames
                rows = list(reader)
            if not fields or "source_object" not in fields or len(rows) != 2:
                raise ValueError("Changed assay manifest structure")
            if origin == run and any(r["source_object"] != "data/retained/cell_target_counts.tsv.gz" for r in rows):
                raise ValueError("Fresh audit must declare its sparse-count source")
            scientific_fields = [field for field in fields if field != "source_object"]
            with (temporary / (name + ".tsv")).open("w", newline="") as handle:
                writer = csv.DictWriter(handle, scientific_fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows({field: row[field] for field in scientific_fields} for row in rows)
        result = reproduction.compare_tsv(temporary / "reference.tsv", temporary / "fresh.tsv")
    return {**result, "explicitly_different_provenance_field": "source_object",
            "scope": "All other fields, including full-assay feature counts and target counts, compared"}


def verify(root: Path, run: Path) -> dict:
    root, run = root.resolve(), run.resolve()
    if run == root or not run.is_dir():
        raise ValueError("Supply a separate completed cell-audit run directory")
    with (run / "cell_audit_replay.json").open() as handle:
        record = json.load(handle)
    if record["rebuild_script_sha256"] != reproduction.digest(root / "scripts/rebuild_cell_audits.R"):
        raise ValueError("Fresh run belongs to a different rebuilder version")
    metadata = "results/tables/annotation/annotated_cell_metadata.tsv.gz"
    if (reproduction.digest(run / metadata) != reproduction.digest(root / metadata) or
            record["metadata_sha256"] != reproduction.digest(root / metadata)):
        raise ValueError("Fresh metadata must match the retained source exactly")
    for source in record["input_files"]:
        path = root / "data/retained" / source["path"]
        if (Path(source["path"]).name != source["path"] or not path.is_file() or
                path.stat().st_size != source["bytes"] or reproduction.digest(path) != source["sha256"]):
            raise ValueError("Fresh audit input differs from the supplied retained counts")
    expected = dict(n_cells=14846, n_union_targets=165, n_expression_summaries=11376,
                    n_unexpected_cell_records=296920, n_positive_QC_records=88779,
                    n_QC_concentration_groups=700, n_biology_summary_groups=1296,
                    n_biology_cell_records=534456, n_biology_donor_contrasts=648,
                    n_biology_share_groups=1296)
    if any(record.get(key) != value for key, value in expected.items()):
        raise ValueError("Incomplete fresh cell-audit execution scope")
    reference_hashes = {name: reproduction.digest(root / "analysis" / name) for name in SCIENCE_OUTPUTS if name not in DERIVED_CELLS}
    comparisons = []
    for name in SCIENCE_OUTPUTS:
        if name in DERIVED_CELLS:
            continue
        result = reproduction.compare_tsv(root / "analysis" / name, run / "analysis" / name)
        comparisons.append({"file": name, **result})
    manifest = assay_manifest_comparison(root, run)
    sparse_expansions = verify_sparse_expansions(root, run)
    ten_genes = verify_cell_expression.verify(run)
    eighteen_genes = verify_biology_cells.verify(run)
    if any(reproduction.digest(root / "analysis" / name) != sha for name, sha in reference_hashes.items()):
        raise ValueError("A retained reference changed during verification")
    # Publish validated generated expansions for subsequent models/verifiers.
    for name in DERIVED_CELLS:
        shutil.copy2(run / "analysis" / name, root / "analysis" / name)
    return dict(sparse_expansions=sparse_expansions, status="PASS", fresh_run=str(run), full_scientific_output_comparisons=comparisons,
                assay_manifest_scientific_fields=manifest,
                ten_gene_independent_verification=ten_genes,
                eighteen_gene_independent_verification=eighteen_genes,
                relative_tolerance=reproduction.REL_TOL, absolute_tolerance=reproduction.ABS_TOL,
                retained_reference_hashes=reference_hashes,
                scope="All 158-gene summary fields compared with retained reference; four expanded cell tables joined directly to authentic sparse inputs; independent normalization/distribution checks on fixed ten/eighteen-gene records; no new raw single-cell reconstruction or biological inference")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    arguments = parser.parse_args()
    result = verify(arguments.root, arguments.run)
    text = json.dumps(result, indent=2) + "\n"
    if arguments.report:
        arguments.report.parent.mkdir(parents=True, exist_ok=True)
        arguments.report.write_text(text)
    print(text, end="")
