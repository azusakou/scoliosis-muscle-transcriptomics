#!/usr/bin/env python3
from __future__ import annotations

from project_paths import raw_data_root, executable
from input_integrity import (digest, validate_pairs, validate_fastq_pair, validate_quantification,
                             require_provenance, write_provenance, directory_fingerprint,
                             validate_sra_paired_counts, require_complete_fastq_counts)

import csv
import json
import tempfile
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_ROOT = raw_data_root()
EXTERNAL_ROOT = RAW_ROOT.parent
RAW = RAW_ROOT / "bulk_gse311180"
FASTQ = RAW / "fastq"
FASTERQ_TMP = RAW / "fasterq_tmp"
QUANT = EXTERNAL_ROOT / "processed" / "gse311180" / "salmon"
INDEX = EXTERNAL_ROOT / "reference" / "gencode_v49" / "salmon_2.3.4_index"
MANIFEST = ROOT / "metadata/audit/gse311180_requant_manifest.tsv"
DOWNLOAD_LOG = ROOT / "logs/gse311180_sra_download.tsv"
STATUS = ROOT / "logs/gse311180_requant_status.tsv"
FASTERQ = executable("fasterq-dump", "AIS_FASTERQ")
SRA_STAT = executable("sra-stat", "AIS_SRA_STAT")
SALMON = executable("salmon", "AIS_SALMON")
if not RAW_ROOT.exists():
    raise SystemExit(f"Configured AIS_RAW_ROOT does not exist: {RAW_ROOT}")
FASTQ.mkdir(parents=True, exist_ok=True)
FASTERQ_TMP.mkdir(parents=True, exist_ok=True)
QUANT.mkdir(parents=True, exist_ok=True)

with MANIFEST.open() as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))
validate_pairs(manifest, expected_pairs=7)
if len({row["sra_run"] for row in manifest}) != len(manifest):
    raise SystemExit("Duplicate SRA accession in requantification manifest")
STATUS.parent.mkdir(parents=True, exist_ok=True)
with DOWNLOAD_LOG.open() as handle:
    downloaded = {r["sra_run"]: r for r in csv.DictReader(handle, delimiter="\t")}
if not all(downloaded.get(r["sra_run"], {}).get("status") == "PASS" for r in manifest):
    raise SystemExit("SRA download/checksum gate is not PASS for all 14 samples")
version = subprocess.run([SALMON, '--version'], check=True, capture_output=True, text=True).stdout.strip()
if version != 'salmon 2.3.4':
    raise SystemExit(f'This recipe requires Salmon 2.3.4; found {version}')
index_stamp = INDEX.with_name(INDEX.name + '.provenance.json')
if not index_stamp.is_file():
    raise SystemExit('Unverified/missing Salmon index; run scripts/07_prepare_gencode_salmon.py')
index_metadata = json.loads(index_stamp.read_text())
if index_metadata['outputs'] != directory_fingerprint(INDEX):
    raise SystemExit('Salmon index content checksum mismatch')
index_sha256 = digest(index_stamp)
# Recheck actual SRA objects and spot totals instead of trusting a historic PASS log.
sra_spots = {}
for row in manifest:
    entry = downloaded[row['sra_run']]
    path = Path(entry['path'])
    if path.stat().st_size != int(row['sra_size_bytes']) or digest(path, 'md5') != row['sra_md5']:
        raise SystemExit(f'SRA size/MD5 mismatch: {row["sra_run"]}')
    summary = subprocess.check_output([SRA_STAT, '--quick', '--xml', str(path)], text=True)
    sra_spots[row['sample_id']] = validate_sra_paired_counts(summary, int(row['read_length']))['spot_count']
fasterq_version = subprocess.run([FASTERQ, '--version'], check=True, capture_output=True, text=True).stdout.strip()
fastq_outputs = {}


def extract(row: dict[str, str]) -> tuple[str, str]:
    sample = row["sample_id"]
    r1, r2 = FASTQ / f"{sample}_1.fastq", FASTQ / f"{sample}_2.fastq"
    stamp = FASTQ / f'{sample}.provenance.json'
    inputs = {'sample_id': sample, 'sra_run': row['sra_run'], 'sra_md5': row['sra_md5'],
              'fasterq_version': fasterq_version, 'split_files': True}
    if r1.exists() or r2.exists():
        saved = require_provenance(stamp, inputs)
        outputs = validate_fastq_pair(r1, r2)
        require_complete_fastq_counts(outputs, sra_spots[sample])
        if saved['outputs'] != outputs:
            raise ValueError(f'FASTQ output checksum mismatch: {sample}')
        fastq_outputs[sample] = outputs
        return sample, "FASTQ_PRESENT"
    sra = downloaded[row["sra_run"]]["path"]
    sample_tmp = FASTERQ_TMP / sample
    sample_tmp.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=FASTQ, prefix=f'{sample}_extract_') as stage:
        proc = subprocess.run([
            FASTERQ, sra, "--split-files", "--threads", "6", "--outdir", stage,
            "--outfile", sample, "--temp", str(sample_tmp)
        ], text=True, capture_output=True)
        if proc.returncode != 0:
            return sample, "FASTERQ_FAILED:" + (proc.stderr or proc.stdout)[-500:]
        staged1, staged2 = Path(stage) / r1.name, Path(stage) / r2.name
        outputs = validate_fastq_pair(staged1, staged2)
        require_complete_fastq_counts(outputs, sra_spots[sample])
        staged1.replace(r1)
        staged2.replace(r2)
    write_provenance(stamp, inputs, outputs)
    fastq_outputs[sample] = outputs
    return sample, "FASTQ_PASS"


extraction: dict[str, str] = {}
# fasterq performs heavy read/write traffic on the same removable volume. Use
# one extraction at a time to avoid USB/exFAT I/O saturation and preserve a
# restartable, checksum-gated workflow.
with ThreadPoolExecutor(max_workers=1) as pool:
    futures = {pool.submit(extract, row): row["sample_id"] for row in manifest}
    for future in as_completed(futures):
        sample, status = future.result()
        extraction[sample] = status
        print(f"{sample}\t{status}", flush=True)
if not all(v in {"FASTQ_PRESENT", "FASTQ_PASS"} for v in extraction.values()):
    raise SystemExit("At least one fasterq-dump extraction failed")


def quantify(row: dict[str, str]) -> tuple[str, str, str]:
    sample = row["sample_id"]
    out = QUANT / sample
    quant = out / "quant.sf"
    stamp = QUANT / f'{sample}.provenance.json'
    inputs = {'sample_id': sample, 'sra_run': row['sra_run'], 'fastq': fastq_outputs[sample],
              'index_provenance_sha256': index_sha256, 'salmon_version': version,
              'libType': 'A', 'validateMappings': True, 'gcBias': True}
    if out.exists():
        saved = require_provenance(stamp, inputs)
        outputs = validate_quantification(out, expected_read_pairs=fastq_outputs[sample]['read_pairs'])
        if saved['outputs'] != outputs:
            raise ValueError(f'Salmon cache checksum mismatch: {sample}')
        return sample, "QUANT_PRESENT", str(quant)
    with tempfile.TemporaryDirectory(dir=QUANT, prefix=f'{sample}_quant_') as stage:
        staged_out = Path(stage) / 'quant'
        proc = subprocess.run([
            SALMON, "quant", "--index", str(INDEX), "--libType", "A",
            "--mates1", str(FASTQ / f"{sample}_1.fastq"),
            "--mates2", str(FASTQ / f"{sample}_2.fastq"),
            "--output", str(staged_out), "--threads", "6", "--validateMappings", "--gcBias"
        ], text=True, capture_output=True)
        if proc.returncode != 0:
            return sample, "SALMON_FAILED:" + (proc.stderr or proc.stdout)[-500:], ""
        outputs = validate_quantification(staged_out, expected_read_pairs=fastq_outputs[sample]['read_pairs'])
        staged_out.replace(out)
    write_provenance(stamp, inputs, outputs)
    return sample, "QUANT_PASS", str(quant)


quantification: dict[str, tuple[str, str]] = {}
with ThreadPoolExecutor(max_workers=2) as pool:
    futures = {pool.submit(quantify, row): row["sample_id"] for row in manifest}
    for future in as_completed(futures):
        sample, status, path = future.result()
        quantification[sample] = (status, path)
        print(f"{sample}\t{status}", flush=True)

rows = []
for row in manifest:
    sample = row["sample_id"]
    qstatus, qpath = quantification.get(sample, ("NOT_RUN", ""))
    rows.append({
        "sample_id": sample, "patient_id": row["patient_id"], "side": row["side"],
        "sra_run": row["sra_run"], "fastq_status": extraction.get(sample, "NOT_RUN"),
        "quant_status": qstatus, "quant_sf": qpath,
    })
with STATUS.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
if not all(r["quant_status"] in {"QUANT_PRESENT", "QUANT_PASS"} for r in rows):
    raise SystemExit(f"Salmon quantification incomplete; see {STATUS}")
print("All 14 GSE311180 paired samples passed FASTQ extraction and Salmon quantification")
