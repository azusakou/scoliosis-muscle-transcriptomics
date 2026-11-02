#!/usr/bin/env python3
from __future__ import annotations

from project_paths import raw_data_root, executable
from input_integrity import digest, validate_download, require_provenance, write_provenance, directory_fingerprint

import csv
import gzip
import hashlib
import os
import re
import subprocess
import tempfile
import urllib.request
from pathlib import Path

RAW_ROOT = raw_data_root()
REF = RAW_ROOT.parent / "reference" / "gencode_v49"
FASTA = REF / "gencode.v49.transcripts.fa.gz"
GTF = REF / "gencode.v49.annotation.gtf.gz"
TX2GENE = REF / "gencode.v49.tx2gene.tsv"
INDEX = REF / "salmon_2.3.4_index"
SALMON = executable("salmon", "AIS_SALMON")
ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "metadata/audit/gencode_v49_reference.tsv"

REF.mkdir(parents=True, exist_ok=True)
EXPECTED = {
    FASTA.name: '08ee7ee86514061de6fd62dc643c2fe962fbc393080589ebc491bd11e22b119e',
    GTF.name: 'd6e6fe0515c95b2a8cd36a853c1989cee9115c736c60237c56ae92b9daaaf7c4',
}
version = subprocess.run([SALMON, '--version'], check=True, capture_output=True, text=True).stdout.strip()
if version != 'salmon 2.3.4':
    raise SystemExit(f'This reference recipe requires Salmon 2.3.4; found {version}')


def download_if_missing(url: str, path: Path) -> None:
    if path.exists():
        validate_download(path)
        if digest(path) != EXPECTED[path.name]:
            raise ValueError(f'GENCODE v49 reference checksum mismatch: {path}')
        return
    request = urllib.request.Request(url, headers={"User-Agent": "ais-paraspinal-ecology/1.0"})
    temporary = None
    try:
        with urllib.request.urlopen(request, timeout=180) as response:
            with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
                temporary = Path(handle.name)
                while True:
                    chunk = response.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
        validate_download(temporary, path.name)
        if digest(temporary) != EXPECTED[path.name]:
            raise ValueError(f'GENCODE v49 reference checksum mismatch: {path}')
        temporary.replace(path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


download_if_missing(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz",
    FASTA,
)
download_if_missing(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz",
    GTF,
)


attribute = re.compile(r'(\w+) "([^"]+)";')
mapping: dict[str, tuple[str, str]] = {}
with gzip.open(GTF, "rt") as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[2] != "transcript":
            continue
        attrs = dict(attribute.findall(fields[8]))
        transcript = attrs.get("transcript_id")
        gene_id = attrs.get("gene_id")
        gene_name = attrs.get("gene_name")
        if transcript and gene_id and gene_name:
            if transcript in mapping and mapping[transcript] != (gene_id, gene_name):
                raise ValueError(f'Conflicting transcript mapping: {transcript}')
            mapping[transcript] = (gene_id, gene_name)

with TX2GENE.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["transcript_id", "gene_id", "gene_name"])
    for transcript, (gene_id, gene_name) in sorted(mapping.items()):
        writer.writerow([transcript, gene_id, gene_name])

index_inputs = {'reference_sha256': EXPECTED, 'salmon_version': version,
                'gencode': True, 'kmer_length': 31, 'mode': 'transcriptome_only'}
index_stamp = REF / 'salmon_2.3.4_index.provenance.json'
if INDEX.exists():
    saved = require_provenance(index_stamp, index_inputs)
    if saved['outputs'] != directory_fingerprint(INDEX):
        raise ValueError('Salmon index checksum mismatch; preserve/move the directory and rebuild')
else:
    with tempfile.TemporaryDirectory(dir=REF, prefix='index_build_') as stage:
        staged_index = Path(stage) / 'index'
        subprocess.run([
            SALMON, "index", "--transcripts", str(FASTA), "--index", str(staged_index),
            "--gencode", "--kmerLen", "31", "--threads", "8"
        ], check=True)
        for name in ('info.json', 'index.ssi', 'refseq.bin'):
            if not (staged_index / name).is_file() or (staged_index / name).stat().st_size == 0:
                raise ValueError(f'Incomplete Salmon index: {name}')
        staged_index.replace(INDEX)
    write_provenance(index_stamp, index_inputs, directory_fingerprint(INDEX))

AUDIT.parent.mkdir(parents=True, exist_ok=True)
with AUDIT.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["field", "value"])
    writer.writerows([
        ("reference", "GENCODE_human_release_49"),
        ("transcript_fasta", str(FASTA)),
        ("annotation_gtf", str(GTF)),
        ("transcript_fasta_sha256", digest(FASTA)),
        ("annotation_gtf_sha256", digest(GTF)),
        ("tx2gene_sha256", digest(TX2GENE)),
        ("n_transcripts_in_tx2gene", len(mapping)),
        ("salmon_version", version),
        ("index_mode", "transcriptome_only_gencode_k31"),
    ])
print(f"GENCODE v49 and Salmon index ready: {len(mapping)} transcripts")
