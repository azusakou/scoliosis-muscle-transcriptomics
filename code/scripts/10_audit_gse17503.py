#!/usr/bin/env python3
"""Recover and gate the five GSE17503 biopsy/culture patient pairs."""

from __future__ import annotations

from project_paths import raw_data_root

import csv
import gzip
import os
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = raw_data_root()
SOFT = RAW_ROOT / "reference_gse17503" / "GSE17503_family.soft.gz"
OUT = ROOT / "metadata/audit/gse17503_samples.tsv"
STATUS = ROOT / "metadata/audit/gse17503_audit_status.tsv"
if not SOFT.exists():
    raise SystemExit(f"Missing {SOFT}")

samples: list[dict[str, str]] = []
current: dict[str, str] | None = None
with gzip.open(SOFT, "rt", encoding="utf-8", errors="replace") as handle:
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE = "):
            if current:
                samples.append(current)
            current = {"sample_id": line.split(" = ", 1)[1]}
        elif current is not None and line.startswith("!Sample_title = "):
            current["title"] = line.split(" = ", 1)[1]
        elif current is not None and line.startswith("!Sample_source_name_ch1 = "):
            current["source_name"] = line.split(" = ", 1)[1]
if current:
    samples.append(current)

rows = []
for sample in samples:
    title = sample.get("title", "")
    match = re.fullmatch(r"(B\d+)\s+[Ii]n\s+(Vitro|Vivo)", title)
    if not match:
        continue
    rows.append({
        "sample_id": sample["sample_id"],
        "title": title,
        "patient_id": match.group(1),
        "material": "culture" if match.group(2).lower() == "vitro" else "biopsy",
        "source_name": sample.get("source_name", ""),
    })

patients = sorted({row["patient_id"] for row in rows})
complete = [
    patient for patient in patients
    if {row["material"] for row in rows if row["patient_id"] == patient} == {"culture", "biopsy"}
]
gate = len(rows) == 10 and len(complete) == 5
if not gate:
    raise SystemExit(f"GSE17503 gate failed: {len(rows)} samples, {len(complete)} complete pairs")

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
with STATUS.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["gate", "status", "value", "interpretation"], delimiter="\t")
    writer.writeheader()
    writer.writerow({"gate": "biopsy_culture_pairing", "status": "PASS", "value": "5_complete_pairs",
                     "interpretation": "Eligible for tissue-versus-culture reference sensitivity only"})
    writer.writerow({"gate": "side_labels", "status": "FAIL", "value": "absent",
                     "interpretation": "Never eligible for concave-versus-convex inference"})
print("GSE17503 audit PASS: five biopsy/culture pairs; side-validation gate FAIL")
