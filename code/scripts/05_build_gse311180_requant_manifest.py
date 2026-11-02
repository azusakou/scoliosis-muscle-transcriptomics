#!/usr/bin/env python3
from __future__ import annotations

from project_paths import raw_data_root
from input_integrity import validate_pairs

import csv
import os
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = raw_data_root()
RAW = RAW_ROOT / "bulk_gse311180"
XML = RAW / "PRJNA1368770_sra_metadata.xml"
AUDIT = ROOT / "metadata/audit/gse311180_samples.tsv"
OUT = ROOT / "metadata/audit/gse311180_requant_manifest.tsv"

with AUDIT.open() as handle:
    audited = list(csv.DictReader(handle, delimiter="\t"))
eligible = {
    row["geo_accession"]: row
    for row in audited
    if row["analysis_role"] == "raw_FASTQ_requantification_candidate"
}

rows: list[dict[str, str | int]] = []
tree = ET.parse(XML)
for package in tree.getroot().findall("EXPERIMENT_PACKAGE"):
    sample = package.find("SAMPLE")
    if sample is None:
        continue
    gsm = ""
    for ext in sample.findall("./IDENTIFIERS/EXTERNAL_ID"):
        if ext.attrib.get("namespace") == "GEO":
            gsm = ext.text or ""
    if gsm not in eligible:
        continue
    runs = package.findall("./RUN_SET/RUN")
    if len(runs) != 1:
        raise SystemExit(f"Expected one run for {gsm}, observed {len(runs)}")
    layout = package.find('./EXPERIMENT/DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/PAIRED')
    if layout is None:
        raise SystemExit(f'Expected paired-read layout for {gsm}')
    run = runs[0]
    sra_files = [f for f in run.findall("./SRAFiles/SRAFile")
                 if f.attrib.get("semantic_name") == "SRA Normalized"]
    if len(sra_files) != 1:
        raise SystemExit(f"Expected one normalized SRA object for {gsm}")
    a = eligible[gsm]
    rows.append({
        "sample_id": a["sample_id"],
        "patient_id": a["patient_id"],
        "side": a["side"],
        "geo_accession": gsm,
        "sra_run": run.attrib["accession"],
        "sra_url": sra_files[0].attrib["url"],
        "sra_size_bytes": int(sra_files[0].attrib["size"]),
        "sra_md5": sra_files[0].attrib["md5"],
        "layout": "PAIRED",
        "read_length": 150,
    })

rows.sort(key=lambda r: (str(r["patient_id"]), str(r["side"])))
if len(rows) != 14 or len({r["patient_id"] for r in rows}) != 7:
    raise SystemExit(f"Expected 14 samples from 7 pairs; observed {len(rows)} samples")
validate_pairs(rows, expected_pairs=7)
if len({row['sra_run'] for row in rows}) != len(rows):
    raise SystemExit('Duplicate SRA accession in manifest')
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

total_gib = sum(int(r["sra_size_bytes"]) for r in rows) / 1024**3
print(f"Wrote {OUT}: 14 samples, 7 pairs, {total_gib:.1f} GiB normalized SRA")
