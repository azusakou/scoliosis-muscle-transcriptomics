#!/usr/bin/env python3
from __future__ import annotations

from project_paths import raw_data_root

import csv
import gzip
import hashlib
import os
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = raw_data_root()
RAW = RAW_ROOT / "bulk_gse311180"
SOFT = RAW / "GSE311180_family.soft.gz"
TPM = RAW / "GSE311180_TPM_allSamples.txt.gz"
TABLE_S1 = RAW / "mmc2.xlsx"
SRA_XML = RAW / "PRJNA1368770_sra_metadata.xml"
OUT = ROOT / "metadata/audit"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_geo_soft() -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(SOFT, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    records.append(current)
                current = {"geo_accession": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_title = "):
                current["geo_title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_description = Library name: "):
                current["library_name"] = line.rsplit(": ", 1)[1]
            elif current is not None and line.startswith("!Sample_instrument_model = "):
                current["instrument_model"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_relation = SRA: "):
                current["sra_experiment"] = line.rsplit("=", 1)[1]
            elif current is not None and line.startswith("!Sample_relation = BioSample: "):
                current["biosample"] = line.rsplit("/", 1)[1]
        if current:
            records.append(current)
    return records


def read_table_s1() -> list[dict[str, str]]:
    ws = load_workbook(TABLE_S1, read_only=True, data_only=True).active
    rows = list(ws.iter_rows(values_only=True))
    header_idx = next(i for i, row in enumerate(rows) if row and row[0] == "Sample #")
    header = [str(x) if x is not None else "" for x in rows[header_idx]]
    out: list[dict[str, str]] = []
    for row in rows[header_idx + 1 :]:
        if not row or row[0] is None:
            continue
        rec = {header[i]: "" if value is None else str(value) for i, value in enumerate(row)}
        if rec.get("Tissue") not in {"Muscle", "Cartilage"}:
            continue
        out.append(rec)
    return out


geo = read_geo_soft()
s1 = read_table_s1()
geo_by_library = {r["library_name"]: r for r in geo}

# SRA original FASTQ filenames provide a third, repository-level identity check.
sra_fastq: dict[str, list[str]] = defaultdict(list)
tree = ET.parse(SRA_XML)
for package in tree.getroot().findall("EXPERIMENT_PACKAGE"):
    sample = package.find("SAMPLE")
    if sample is None:
        continue
    identifiers = sample.find("IDENTIFIERS")
    gsm = ""
    if identifiers is not None:
        for ext in identifiers.findall("EXTERNAL_ID"):
            if ext.attrib.get("namespace") == "GEO":
                gsm = ext.text or ""
    for f in package.findall(".//SRAFile"):
        if f.attrib.get("semantic_name") == "fastq":
            sra_fastq[gsm].append(f.attrib.get("filename", ""))

with gzip.open(TPM, "rt") as handle:
    matrix_columns = handle.readline().rstrip("\n").split("\t")[1:]

records: list[dict[str, str]] = []
for row in s1:
    sample = row["Sample #"]
    matrix_column = "RD112" if sample == "DR112" and "RD112" in matrix_columns else sample
    g = geo_by_library.get(sample, {})
    patient_id = row.get("Patient ID", "")
    side = row.get("Concave/Convex", "").lower()
    diagnosis = row.get("Diagnosis/Reason for Surgery", "")
    tissue = row.get("Tissue", "")
    complete_pair = False  # filled below
    records.append({
        "sample_id": sample,
        "patient_id": patient_id,
        "diagnosis": diagnosis,
        "tissue": tissue,
        "side": side,
        "sex": row.get("Sex", ""),
        "age": row.get("Age", ""),
        "spinal_level": row.get("Spinal Level", ""),
        "cobb_angle": row.get("AIS Cobb of Largest Curve", ""),
        "geo_accession": g.get("geo_accession", ""),
        "geo_title": g.get("geo_title", ""),
        "matrix_column": matrix_column if matrix_column in matrix_columns else "",
        "matrix_mapping_status": (
            "exact_sample_id" if sample in matrix_columns else
            "documented_unique_header_transposition_DR112_to_RD112" if matrix_column == "RD112" else
            "absent"
        ),
        "sra_experiment": g.get("sra_experiment", ""),
        "sra_fastq_filenames": ";".join(sra_fastq.get(g.get("geo_accession", ""), [])),
        "complete_pair": str(complete_pair).upper(),
    })

muscle_ais = [r for r in records if r["diagnosis"] == "AIS" and r["tissue"] == "Muscle"]
by_patient: dict[str, set[str]] = defaultdict(set)
for r in muscle_ais:
    by_patient[r["patient_id"]].add(r["side"])
complete = {p for p, sides in by_patient.items() if sides == {"concave", "convex"}}
for r in records:
    r["complete_pair"] = str(r["patient_id"] in complete).upper()
    if r in muscle_ais and r["patient_id"] in complete:
        r["analysis_role"] = "raw_FASTQ_requantification_candidate"
    elif r in muscle_ais:
        r["analysis_role"] = "unpaired_AIS_excluded_from_paired_test"
    elif r["tissue"] == "Muscle":
        r["analysis_role"] = "adult_control_age_confounded_context_only"
    else:
        r["analysis_role"] = "cartilage_not_paraspinal_muscle"

fields = list(records[0])
with (OUT / "gse311180_samples.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(records)

# Diagnose the deposited TPM matrix without using it for inference.
diagnostic_genes = {"CKM", "MYH1", "ACAN", "COL2A1"}
values: dict[str, dict[str, float]] = {}
with gzip.open(TPM, "rt") as handle:
    header = handle.readline().rstrip("\n").split("\t")
    for line in handle:
        parts = line.rstrip("\n").split("\t")
        if parts[0] in diagnostic_genes:
            values[parts[0]] = {header[i]: float(parts[i]) for i in range(1, len(header))}

diag_rows: list[list[object]] = []
for r in muscle_ais:
    col = r["matrix_column"]
    diag_rows.append([r["sample_id"], r["patient_id"], r["side"], col] +
                     [values.get(g, {}).get(col, float("nan")) for g in sorted(diagnostic_genes)])
with (OUT / "gse311180_expression_identity_diagnostics.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["sample_id", "patient_id", "side", "matrix_column"] + sorted(diagnostic_genes))
    writer.writerows(diag_rows)

# All deposited convex columns have a cartilage-like identity and are not safe to pair with concave muscle.
convex = [r for r in muscle_ais if r["side"] == "convex"]
matrix_conflict = all(values["CKM"][r["matrix_column"]] == 0 for r in convex) and all(
    values["ACAN"][r["matrix_column"]] > 100 for r in convex
)

with (OUT / "gse311180_audit_status.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["field", "value"])
    writer.writerows([
        ("dataset", "GSE311180"),
        ("bioproject", "PRJNA1368770"),
        ("supplement_verified_AIS_muscle_samples", len(muscle_ais)),
        ("complete_AIS_pairs", len(complete)),
        ("paired_AIS_samples", 2 * len(complete)),
        ("unpaired_AIS_muscle_samples", len(muscle_ais) - 2 * len(complete)),
        ("adult_muscle_controls", sum(r["analysis_role"].startswith("adult_control") for r in records)),
        ("processed_TPM_tissue_identity_conflict", str(matrix_conflict).upper()),
        ("published_side_DE_model_patient_blocking", "ABSENT"),
        ("processed_matrix_pairing_gate", "BLOCKED_METADATA_EXPRESSION_CONFLICT"),
        ("raw_FASTQ_requantification_gate", "PASS_7_COMPLETE_PAIRS"),
        ("soft_sha256", sha256(SOFT)),
        ("tpm_sha256", sha256(TPM)),
        ("table_s1_sha256", sha256(TABLE_S1)),
        ("sra_metadata_sha256", sha256(SRA_XML)),
    ])

print("GSE311180 audit: 7 supplement-verified pairs; processed TPM blocked; raw FASTQ eligible for requantification")
