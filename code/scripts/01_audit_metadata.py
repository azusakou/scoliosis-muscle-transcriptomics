#!/usr/bin/env python3
"""Create donor/side audit tables and enforce the pre-analysis design gate."""

from __future__ import annotations

from project_paths import raw_data_root

import argparse
import csv
import gzip
import hashlib
import json
import pathlib
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict


DEFAULT_ROOT = raw_data_root()


def text(node, path, default=""):
    found = node.find(path)
    return found.text.strip() if found is not None and found.text else default


def parse_sra(path: pathlib.Path):
    root = ET.parse(path).getroot()
    rows = []
    for package in root.findall("EXPERIMENT_PACKAGE"):
        experiment = package.find("EXPERIMENT")
        sample = package.find("SAMPLE")
        run = package.find("RUN_SET/RUN")
        if experiment is None or sample is None or run is None:
            continue
        attrs = {}
        for item in sample.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
            attrs[text(item, "TAG")] = text(item, "VALUE")
        title = text(experiment, "TITLE")
        description = text(experiment, "DESIGN/DESIGN_DESCRIPTION")
        side = "concave" if description.lower().startswith("concave") else "convex" if description.lower().startswith("convex") else "unknown"
        patient_match = re.match(r"([^_]+)_", title)
        patient_id = patient_match.group(1) if patient_match else "unknown"
        biosample = ""
        for external in sample.findall("IDENTIFIERS/EXTERNAL_ID"):
            if external.attrib.get("namespace") == "BioSample":
                biosample = external.text or ""
        rows.append({
            "dataset": "PRJNA722100",
            "patient_id": patient_id,
            "side": side,
            "experiment_title": title,
            "design_description": description,
            "experiment_accession": experiment.attrib.get("accession", ""),
            "run_accession": run.attrib.get("accession", ""),
            "sra_sample_accession": sample.attrib.get("accession", ""),
            "biosample_accession": biosample,
            "sample_alias": sample.attrib.get("alias", ""),
            "age": attrs.get("age", ""),
            "sex": attrs.get("sex", ""),
            "tissue_submitter": attrs.get("tissue", ""),
            "organism": text(sample, "SAMPLE_NAME/SCIENTIFIC_NAME"),
            "library_strategy": text(experiment, "DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY"),
            "library_source": text(experiment, "DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE"),
            "platform": text(experiment, "PLATFORM/ILLUMINA/INSTRUMENT_MODEL"),
            "total_spots": run.attrib.get("total_spots", ""),
            "total_bases": run.attrib.get("total_bases", ""),
        })
    return sorted(rows, key=lambda row: (row["patient_id"], row["side"]))


def parse_geo_soft(path: pathlib.Path):
    series = defaultdict(list)
    samples = []
    current = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    samples.append(current)
                current = {"soft_sample_id": line.split("=", 1)[1].strip(), "characteristics": [], "relations": []}
                continue
            if current is None:
                if line.startswith("!Series_") and " = " in line:
                    key, value = line.split(" = ", 1)
                    series[key[8:]].append(value)
                continue
            if line.startswith("^SAMPLE"):
                continue
            if line.startswith("!") and " = " in line:
                key, value = line.split(" = ", 1)
                field = key[8:] if key.startswith("!Sample_") else key
                if field == "characteristics_ch1":
                    current["characteristics"].append(value)
                elif field == "relation":
                    current["relations"].append(value)
                elif field not in current:
                    current[field] = value
    if current:
        samples.append(current)

    rows = []
    side_map = {"A": "concave", "T": "convex"}
    for sample in samples:
        title = sample.get("title", "")
        processed_match = re.search(r"\[([^\]]+) in processed data\]", title)
        processed_id = processed_match.group(1) if processed_match else title.split()[0]
        id_match = re.fullmatch(r"(IS|CS)(\d+)([AT])", processed_id)
        status = ""
        for item in sample["characteristics"]:
            if item.lower().startswith("subject status:"):
                status = item.split(":", 1)[1].strip()
        relation_map = {}
        for relation in sample["relations"]:
            if ": " in relation:
                key, value = relation.split(": ", 1)
                relation_map[key] = value.rsplit("/", 1)[-1].split("=")[-1]
        if id_match:
            cohort_code, donor_number, side_code = id_match.groups()
            patient_id = f"{cohort_code}{donor_number}"
            side = side_map[side_code]
        else:
            cohort_code = patient_id = "unknown"
            side = "unknown"
        rows.append({
            "dataset": "GSE254300",
            "geo_accession": sample.get("geo_accession", sample["soft_sample_id"]),
            "submitted_title": title,
            "processed_sample_id": processed_id,
            "cohort": "AIS" if cohort_code == "IS" else "congenital_scoliosis" if cohort_code == "CS" else "unknown",
            "patient_id": patient_id,
            "side": side,
            "subject_status": status,
            "tissue": sample.get("source_name_ch1", ""),
            "instrument_model": sample.get("instrument_model", ""),
            "library_strategy": sample.get("library_strategy", ""),
            "library_source": sample.get("library_source", ""),
            "sra_experiment": relation_map.get("SRA", ""),
            "biosample_accession": relation_map.get("BioSample", ""),
        })
    return series, sorted(rows, key=lambda row: (row["cohort"], row["patient_id"], row["side"]))


def write_tsv(path: pathlib.Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def paired_complete(rows, cohort=None):
    observed = defaultdict(set)
    for row in rows:
        if cohort is None or row.get("cohort") == cohort:
            observed[row["patient_id"]].add(row["side"])
    return bool(observed) and all(sides == {"concave", "convex"} for sides in observed.values()), observed


def xlsx_header(path: pathlib.Path):
    with zipfile.ZipFile(path) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
            for item in root.findall("x:si", namespace):
                shared.append("".join(node.text or "" for node in item.findall(".//x:t", namespace)))
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        namespace = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        first_row = sheet.find(".//x:sheetData/x:row", namespace)
        values = []
        for cell in first_row.findall("x:c", namespace):
            value = text(cell, "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v")
            values.append(shared[int(value)] if cell.attrib.get("t") == "s" else value)
        return values


def matrix_market_header(path: pathlib.Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        metadata = ""
        for line in handle:
            if line.startswith("%metadata_json:"):
                metadata = line.split(":", 1)[1].strip()
            if not line.startswith("%"):
                genes, barcodes, nonzero = map(int, line.split())
                return genes, barcodes, nonzero, metadata
    raise ValueError(f"No MatrixMarket dimensions in {path}")


def sha256(path: pathlib.Path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=pathlib.Path, default=DEFAULT_ROOT)
    parser.add_argument("--output", type=pathlib.Path, default=pathlib.Path("metadata/audit"))
    args = parser.parse_args()
    raw = args.raw_root
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    scrna = parse_sra(raw / "metadata/PRJNA722100_sra_experiments.xml")
    series, bulk = parse_geo_soft(raw / "metadata/GSE254300_family.soft")
    write_tsv(out / "scrna_samples.tsv", scrna)
    write_tsv(out / "bulk_samples.tsv", bulk)

    sc_paired, sc_pairs = paired_complete(scrna)
    ais_paired, ais_pairs = paired_complete(bulk, "AIS")
    cs_paired, cs_pairs = paired_complete(bulk, "congenital_scoliosis")
    count_header = xlsx_header(raw / "bulk/GSE254300_readcount.xlsx")
    expected_bulk_header = ["GeneID"] + [row["processed_sample_id"] for row in bulk if row["cohort"] == "congenital_scoliosis"] + [row["processed_sample_id"] for row in bulk if row["cohort"] == "AIS"]
    header_matches_metadata = set(count_header[1:]) == set(expected_bulk_header[1:])

    matrix_rows = []
    for folder in sorted((raw / "scrna/matrices").glob("*_raw_feature_bc_matrix")):
        matrix = folder / "matrix.mtx.gz"
        if matrix.exists():
            genes, barcodes, nonzero, metadata_json = matrix_market_header(matrix)
            name_match = re.match(r"(.+)_(concave|convex)_raw_feature_bc_matrix", folder.name)
            matrix_rows.append({
                "patient_id": name_match.group(1),
                "side": name_match.group(2),
                "features": genes,
                "raw_barcodes": barcodes,
                "nonzero_entries": nonzero,
                "matrix_metadata": metadata_json,
                "path": str(matrix),
            })
    write_tsv(out / "scrna_matrix_manifest.tsv", matrix_rows)

    feature_counts = sorted({row["features"] for row in matrix_rows})
    gates = {
        "metadata_retrieved": len(scrna) == 4 and len(bulk) == 20,
        "scrna_patient_id_confirmed": all(row["patient_id"] != "unknown" for row in scrna),
        "scrna_side_confirmed": all(row["side"] in {"concave", "convex"} for row in scrna),
        "scrna_pairing_complete": sc_paired and len(sc_pairs) == 2,
        "bulk_patient_id_confirmed": all(row["patient_id"] != "unknown" for row in bulk),
        "bulk_side_confirmed": all(row["side"] in {"concave", "convex"} for row in bulk),
        "bulk_ais_pairing_complete": ais_paired and len(ais_pairs) == 5,
        "bulk_cs_pairing_complete": cs_paired and len(cs_pairs) == 5,
        "bulk_count_header_matches_soft_metadata": header_matches_metadata,
        "scrna_raw_matrices_present": len(matrix_rows) == 4,
        "scrna_feature_reference_consistent": len(feature_counts) == 1,
    }
    status = {
        "datasets": {
            "PRJNA722100": {"samples": len(scrna), "patients": len(sc_pairs), "paired": sc_paired, "platforms": sorted({r["platform"] for r in scrna})},
            "GSE254300": {"samples": len(bulk), "ais_patients": len(ais_pairs), "congenital_patients": len(cs_pairs), "paired": ais_paired and cs_paired, "platforms": sorted({r["instrument_model"] for r in bulk})},
        },
        "gates": gates,
        "analysis_policy": {
            "paired_bulk_inference": "permitted" if gates["bulk_ais_pairing_complete"] else "blocked",
            "paired_scrna_pseudobulk": "exploratory_effect_sizes_only",
            "scrna_reason": "Only two paired AIS donors are present; inferential p-values, stable latent-factor estimation, and majority-donor communication filters are not defensible.",
            "cross_cell_axis": "exploratory and must be validated in independent bulk data",
            "reference_harmonization": "required" if len(feature_counts) > 1 else "not_required",
        },
        "series_overall_design": series.get("overall_design", []),
    }
    (out / "audit_status.json").write_text(json.dumps(status, indent=2) + "\n", encoding="utf-8")

    checksum_targets = []
    for subdir in ("metadata", "bulk", "scrna"):
        for path in sorted((raw / subdir).glob("*")):
            if path.is_file() and not path.name.startswith("._"):
                checksum_targets.append(path)
    with (out / "checksums.sha256").open("w", encoding="utf-8") as handle:
        for path in checksum_targets:
            handle.write(f"{sha256(path)}  {path.relative_to(raw)}\n")

    print(json.dumps(status, indent=2))
    return 0 if all(value for key, value in gates.items() if key != "scrna_feature_reference_consistent") else 2


if __name__ == "__main__":
    sys.exit(main())

