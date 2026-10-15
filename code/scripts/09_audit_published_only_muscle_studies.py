#!/usr/bin/env python3
"""Audit AIS muscle studies that cannot enter patient-level meta-analysis."""

from __future__ import annotations

from project_paths import raw_data_root

import csv
import hashlib
import os
import re
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = raw_data_root()
OUT = ROOT / "metadata/audit/published_only_muscle_studies.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def docx_text(path: Path) -> str:
    if not zipfile.is_zipfile(path):
        raise RuntimeError(f"Not a valid DOCX ZIP package: {path}")
    with zipfile.ZipFile(path) as archive:
        xml = archive.read("word/document.xml").decode("utf-8")
    return " ".join(re.sub(r"<[^>]+>", " ", xml).split())


h19 = RAW / "literature_h19_2018/10020_2018_49_MOESM1_ESM.docx"
tent5a_supp = RAW / "literature_tent5a_2022/CPR-55-e13183-s002.docx"
tent5a_xml = RAW / "literature_tent5a_2022/PMC8891553_fullTextXML.xml"
for required in (h19, tent5a_supp, tent5a_xml):
    if not required.exists() or required.stat().st_size == 0:
        raise SystemExit(f"Missing published-study audit source: {required}")

h19_text = docx_text(h19)
tent5a_supp_text = docx_text(tent5a_supp)
tent5a_article = tent5a_xml.read_text(encoding="utf-8")
if "Table S3. Differentially expressed genes" not in h19_text:
    raise SystemExit("H19 supplement no longer contains the expected DEG-table heading")
if "Twenty female AIS patients" not in tent5a_article:
    raise SystemExit("TENT5A full text no longer supports the audited collection cohort")
if "Supplement Table 1" not in tent5a_supp_text:
    raise SystemExit("TENT5A supplement structure is not as expected")

rows = [
    {
        "study": "H19_ADIPOQ_2018",
        "pmid": "30241458",
        "tissue": "bilateral apical paravertebral muscle",
        "reported_rnaseq_pairs": "5",
        "reported_larger_validation_pairs": "60",
        "published_degs": "40",
        "patient_level_matrix_public": "NO",
        "meta_analysis_eligible": "NO",
        "reason": "No public patient-level matrix; DEG table is an embedded EMF drawing",
        "source_file": str(h19),
        "source_sha256": sha256(h19),
    },
    {
        "study": "TENT5A_2022",
        "pmid": "35137485",
        "tissue": "bilateral apical multifidus",
        "reported_rnaseq_pairs": "NA_NOT_IDENTIFIED",
        "reported_larger_validation_pairs": "20_collected_not_explicitly_all_rnaseq",
        "published_degs": "58",
        "patient_level_matrix_public": "NO",
        "meta_analysis_eligible": "NO",
        "reason": "RNA-seq subset/sample membership and patient-level counts are not public",
        "source_file": str(tent5a_xml),
        "source_sha256": sha256(tent5a_xml),
    },
]
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote {OUT} ({len(rows)} published-only studies)")
