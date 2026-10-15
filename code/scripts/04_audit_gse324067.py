#!/usr/bin/env python3
from __future__ import annotations

from project_paths import raw_data_root

import csv
import gzip
import hashlib
import os
import re
import tarfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = raw_data_root()
RAW = RAW_ROOT / "bulk_gse324067"
SOFT = RAW / "GSE324067_family.soft.gz"
ARCHIVE = RAW / "GSE324067_RAW.tar"
SRA_XML = RAW / "PRJNA1433384_sra_metadata.xml"
OUT = ROOT / "metadata/audit"
OUT.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def side_from_text(text: str) -> str:
    low = text.lower().replace("-", "_")
    if "non_curve" in low or "non curve" in low or "noncurve" in low:
        return "non_curve"
    if "concave" in low or re.search(r"(^|[_-])cc([_-]|$)", low):
        return "concave"
    if "convex" in low or re.search(r"(^|[_-])cv([_-]|$)", low):
        return "convex"
    return "unknown"


# GEO records, retained specifically to expose attachment-label conflicts.
geo: list[dict[str, str]] = []
current: dict[str, str] | None = None
with gzip.open(SOFT, "rt", errors="replace") as handle:
    for raw_line in handle:
        line = raw_line.rstrip("\n")
        if line.startswith("^SAMPLE = "):
            if current:
                geo.append(current)
            current = {"geo_accession": line.split(" = ", 1)[1]}
        elif current is not None and line.startswith("!Sample_title = "):
            current["geo_title"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_description = Library name: "):
            current["library_name"] = line.rsplit(": ", 1)[1]
        elif current is not None and line.startswith("!Sample_supplementary_file_1 = "):
            current["processed_file"] = line.rsplit("/", 1)[1]
    if current:
        geo.append(current)
geo_by_file = {r["processed_file"]: r for r in geo}

# SRA original filenames independently verify the GEO title/library-name side for each GSM.
raw_fastq_sides: dict[str, set[str]] = defaultdict(set)
tree = ET.parse(SRA_XML)
for package in tree.getroot().findall("EXPERIMENT_PACKAGE"):
    sample = package.find("SAMPLE")
    if sample is None:
        continue
    gsm = ""
    for ext in sample.findall("./IDENTIFIERS/EXTERNAL_ID"):
        if ext.attrib.get("namespace") == "GEO":
            gsm = ext.text or ""
    for f in package.findall(".//SRAFile"):
        if f.attrib.get("semantic_name") == "fastq":
            raw_fastq_sides[gsm].add(side_from_text(f.attrib.get("filename", "")))

rows: list[dict[str, str]] = []
with tarfile.open(ARCHIVE) as archive:
    for member in archive.getmembers():
        if not member.name.endswith(".tabular.gz"):
            continue
        extracted = archive.extractfile(member)
        if extracted is None:
            raise RuntimeError(f"Cannot read {member.name}")
        with gzip.open(extracted, "rt") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            first = handle.readline().rstrip("\n").split("\t")
        if len(header) != 2 or len(first) != 2:
            raise RuntimeError(f"Unexpected count format in {member.name}")
        internal = header[1]
        m = re.search(r"SC0*(\d+)", internal, flags=re.I)
        if not m:
            raise RuntimeError(f"No patient ID in {internal}")
        patient = f"SC{int(m.group(1)):02d}"
        processed_side = side_from_text(internal)
        filename_side = side_from_text(member.name)
        g = geo_by_file.get(member.name, {})
        gsm = g.get("geo_accession", "")
        geo_side = side_from_text(g.get("geo_title", ""))
        fastq = raw_fastq_sides.get(gsm, set())
        raw_side = next(iter(fastq)) if len(fastq) == 1 else "unknown"
        rows.append({
            "processed_file": member.name,
            "processed_header": internal,
            "patient_id": patient,
            "processed_side": processed_side,
            "filename_side": filename_side,
            "geo_accession_attached": gsm,
            "geo_title": g.get("geo_title", ""),
            "geo_title_side": geo_side,
            "raw_fastq_side_for_attached_GSM": raw_side,
            "processed_label_agreement": str(processed_side == filename_side).upper(),
            "geo_attachment_conflict": str(processed_side != geo_side or processed_side != raw_side).upper(),
            "analysis_role": "adjacent_tissue_specificity_negative_control",
            "inference_label_source": "processed_file_internal_header_plus_basename",
        })

by_patient: dict[str, set[str]] = defaultdict(set)
for row in rows:
    by_patient[row["patient_id"]].add(row["processed_side"])
complete = {p for p, sides in by_patient.items() if sides == {"concave", "convex", "non_curve"}}
if len(rows) != 18 or len(complete) != 6:
    observed = "; ".join(f"{p}={sorted(s)}" for p, s in sorted(by_patient.items()))
    raise SystemExit(f"Expected 18 files and 6 complete triplets; observed {len(rows)} and {len(complete)} ({observed})")
if not all(r["processed_label_agreement"] == "TRUE" for r in rows):
    raise SystemExit("Processed filename and internal count header disagree")

fields = list(rows[0])
with (OUT / "gse324067_samples.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda r: (r["patient_id"], r["processed_side"])))

with (OUT / "gse324067_audit_status.tsv").open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["field", "value"])
    writer.writerows([
        ("dataset", "GSE324067"),
        ("bioproject", "PRJNA1433384"),
        ("tissue", "primary_osteoblasts_from_spinal_bone"),
        ("patients", len(complete)),
        ("processed_samples", len(rows)),
        ("complete_concave_convex_noncurve_triplets", len(complete)),
        ("geo_attachment_conflicts", sum(r["geo_attachment_conflict"] == "TRUE" for r in rows)),
        ("processed_label_gate", "PASS_INTERNAL_HEADER_AND_BASENAME_AGREE"),
        ("direct_muscle_validation_gate", "FAIL_DIFFERENT_TISSUE"),
        ("analysis_role", "adjacent_tissue_specificity_negative_control"),
        ("soft_sha256", sha256(SOFT)),
        ("raw_archive_sha256", sha256(ARCHIVE)),
        ("sra_metadata_sha256", sha256(SRA_XML)),
    ])

print("GSE324067 audit PASS for adjacent-tissue analysis: 6 complete triplets; GEO attachment conflicts explicitly retained")
