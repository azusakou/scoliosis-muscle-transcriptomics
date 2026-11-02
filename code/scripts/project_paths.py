"""Portable paths and external-tool lookup; importing this module writes nothing."""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Direct read aliases; no files are restored, copied, or overwritten.  Only
# actual producer outputs can create the old paths in a fresh raw workspace.
ANALYSIS_INPUT_ALIASES = {
    "results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz": "analysis/bulk_gate_inputs/original_A_gene_results.tsv.gz",
    "results/tables/bulk_gse311180_requant/GSE311180_paired_edgeR_full.tsv.gz": "analysis/bulk_gate_inputs/original_B_gene_results.tsv.gz",
    "results/tables/bulk_gse311180_requant/GSE311180_QC.tsv": "analysis/bulk_gate_inputs/original_B_QC.tsv",
    "analysis/bulk_edgeR_null_logFC_GSE254300.tsv.gz": "analysis/bulk_gate_inputs/A32_logFC.tsv.gz",
    "analysis/bulk_edgeR_null_grid_GSE254300.tsv": "analysis/bulk_gate_inputs/A32_grid.tsv",
    "analysis/audit_gse311180_rna_expected_pairs.tsv": "analysis/bulk_gate_inputs/RNA_expected_pairs.tsv",
    "analysis/bulk_gate_normalization_factors.tsv": "analysis/bulk_gate_inputs/published_retained_TMM_factors.tsv",
    "results/figures/source_data/Figure_S3_osteoblast_control.tsv": "results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_paired_effects.tsv",
}


def resolve_analysis_input(path, root=PROJECT_ROOT) -> Path:
    """Select fresh producer output before fixed baseline; reject unsafe paths.

    An existing fresh file is always returned, even if its bytes differ.  No
    failed hash or changed result is hidden by falling back to the baseline.
    """
    import sys
    name = str(path)
    relative = Path(name)
    if (not name or relative.is_absolute() or "\\" in name or
            any(part in {"", ".", ".."} for part in name.split("/")) or
            re.match(r"^[A-Za-z]:|^~", name)):
        raise ValueError("Unsafe analysis input path")
    root = Path(root).resolve(strict=True)

    def checked(value):
        candidate = root
        for part in Path(value).parts:
            candidate = candidate / part
            if candidate.is_symlink():
                raise ValueError("Symlink analysis input path: " + value)
        if not candidate.exists():
            return None
        resolved = candidate.resolve(strict=True)
        if root not in resolved.parents or not resolved.is_file():
            raise ValueError("Analysis input escapes root or is not a regular file")
        return resolved

    fresh = checked(name)
    canonical = ANALYSIS_INPUT_ALIASES.get(name)
    fixed = checked(canonical) if canonical else None
    if fresh is not None:
        if canonical:
            print("Analysis input scope: fresh producer output " + name +
                  "; immutable baseline " + canonical + " is not substituted.", file=sys.stderr)
        return fresh
    if fixed is not None:
        return fixed
    raise FileNotFoundError("Missing analysis input: " + name)


def raw_data_root() -> Path:
    return Path(os.environ.get("AIS_RAW_ROOT", str(PROJECT_ROOT / "data/raw/public"))).expanduser().resolve()


def executable(name: str, variable: str) -> str:
    configured = os.environ.get(variable, name)
    found = shutil.which(configured)
    if not found:
        raise RuntimeError(f"Required executable {name!r} not found; install it or set {variable}")
    return found


def curl_retry_options(binary: str) -> list[str]:
    """Use the broad retry flag only on curl versions that implement it."""
    result = subprocess.run([binary, "--version"], check=True,
                            capture_output=True, text=True)
    match = re.match(r"curl (\d+)\.(\d+)\.(\d+)(?:\s|[-.])", result.stdout)
    if not match:
        raise RuntimeError("Cannot determine curl version from --version output")
    version = tuple(map(int, match.groups()))
    options = ["--retry", "8"]
    if version >= (7, 71, 0):
        options.append("--retry-all-errors")
    return options
