"""Input format, paired-read and provenance checks; standard library only."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
import tarfile
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

FASTQ_QUALITY_BYTES = bytes(range(33, 127))


def validate_sra_paired_counts(xml: str, expected_read_length: int) -> dict:
    """Read native SRA totals for this fixed-length, biological paired recipe."""
    if type(expected_read_length) is not int or expected_read_length <= 0 or '<!DOCTYPE' in xml.upper():
        raise ValueError('Invalid SRA read-length/XML declaration')
    try:
        root = ET.fromstring(xml)
        counts = {key: int(root.attrib[key]) for key in ('spot_count', 'base_count', 'base_count_bio')}
    except (ET.ParseError, KeyError, ValueError) as error:
        raise ValueError('Malformed SRA count summary') from error
    if (root.tag != 'Run' or counts['spot_count'] <= 0 or
            counts['base_count'] != counts['base_count_bio'] or
            counts['base_count_bio'] != counts['spot_count'] * 2 * expected_read_length):
        raise ValueError('SRA totals do not match the recorded biological paired-read recipe')
    return counts


def require_complete_fastq_counts(outputs: dict, expected_sra_spots: int) -> None:
    if (type(expected_sra_spots) is not int or expected_sra_spots <= 0 or
            type(outputs.get('read_pairs')) is not int or outputs['read_pairs'] != expected_sra_spots):
        raise ValueError('FASTQ records do not cover every original SRA spot')


def digest(path: Path, algorithm: str = 'sha256') -> str:
    value = hashlib.new(algorithm)
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def validate_download(path: Path, declared_name: str | None = None) -> None:
    name = (declared_name or path.name).lower()
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Empty/missing download: {path}')
    with path.open('rb') as handle:
        prefix = handle.read(1024).lstrip().lower()
    if not name.endswith('.html') and (prefix.startswith(b'<!doctype html') or prefix.startswith(b'<html')):
        raise ValueError(f'HTML error/login page received instead of data: {path}')
    if name.endswith(('.xlsx', '.docx', '.zip')):
        with zipfile.ZipFile(path) as handle:
            bad = handle.testzip()
            if bad:
                raise ValueError(f'ZIP CRC failed: {bad}')
            if name.endswith('.xlsx') and 'xl/workbook.xml' not in handle.namelist():
                raise ValueError(f'Not an XLSX workbook: {path}')
            if name.endswith('.docx') and 'word/document.xml' not in handle.namelist():
                raise ValueError(f'Not a DOCX document: {path}')
    elif name.endswith(('.tar', '.tar.gz', '.tgz')):
        # The public cellranger downloads have .tar.gz names but contain ZIPs.
        if zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as handle:
                bad = handle.testzip()
                if bad:
                    raise ValueError(f'ZIP CRC failed: {bad}')
        else:
            with tarfile.open(path) as handle:
                for member in handle:
                    if member.isfile():
                        source = handle.extractfile(member)
                        if source is None:
                            raise ValueError(f'Unreadable TAR member: {member.name}')
                        with source:
                            while source.read(8 * 1024 * 1024):
                                pass
    elif name.endswith('.gz'):
        # Read to EOF: a valid gzip header alone does not establish completeness.
        with gzip.open(path, 'rb') as handle:
            while handle.read(8 * 1024 * 1024):
                pass
    elif name.endswith('.json'):
        with path.open() as handle:
            json.load(handle)
    elif name.endswith('.xml'):
        ET.parse(path)
    elif name.endswith('.soft'):
        with path.open(errors='strict') as handle:
            if not any(line.startswith('^SERIES = GSE') for line in handle):
                raise ValueError(f'Not a GEO SOFT series: {path}')


def validate_pairs(rows: list[dict], sample_column: str = 'sample_id', expected_pairs: int | None = None) -> None:
    samples, pairs = set(), {}
    for row in rows:
        sample, patient, side = (row.get(key, '') for key in (sample_column, 'patient_id', 'side'))
        if not all(isinstance(value, str) and value.strip() for value in (sample, patient, side)):
            raise ValueError('Missing sample, patient or side')
        if sample in samples or side not in {'convex', 'concave'}:
            raise ValueError(f'Duplicate sample or invalid side: {sample}')
        samples.add(sample)
        pairs.setdefault(patient, []).append(side)
    if not pairs or any(sorted(sides) != ['concave', 'convex'] for sides in pairs.values()):
        raise ValueError('Each patient must have exactly one concave and one convex library')
    if expected_pairs is not None and len(pairs) != expected_pairs:
        raise ValueError(f'Expected {expected_pairs} complete pairs, found {len(pairs)}')


def _read_fastq(handle, path: Path):
    header = handle.readline()
    if not header:
        return None
    sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
    if not (header.startswith(b'@') and sequence and plus.startswith(b'+') and quality):
        raise ValueError(f'Truncated/malformed FASTQ record: {path}')
    sequence, quality = sequence.rstrip(b'\r\n'), quality.rstrip(b'\r\n')
    # C-level deletion checks every byte against the same inclusive 33..126
    # range without a Python generator invocation for each base in large runs.
    if not sequence or len(sequence) != len(quality) or quality.translate(None, FASTQ_QUALITY_BYTES):
        raise ValueError(f'Invalid FASTQ sequence/quality lengths: {path}')
    identifier = header[1:].split()[0]
    if identifier.endswith((b'/1', b'/2')):
        identifier = identifier[:-2]
    return identifier


def validate_fastq_pair(first: Path, second: Path) -> dict:
    n = 0
    # Alternate paired records from large buffers, not alternating 8-KiB reads
    # from two huge files on rotational storage. Validation still reads every
    # record and computes the complete per-mate checksums below.
    with first.open('rb', buffering=16 * 1024 * 1024) as left, second.open('rb', buffering=16 * 1024 * 1024) as right:
        while True:
            a, b = _read_fastq(left, first), _read_fastq(right, second)
            if a is None and b is None:
                break
            if a is None or b is None or a != b:
                raise ValueError(f'Unpaired FASTQ identifiers or unequal record counts: {first}, {second}')
            n += 1
    if n == 0:
        raise ValueError('Empty paired FASTQs')
    return {'read_pairs': n, 'mate1_sha256': digest(first), 'mate2_sha256': digest(second)}


def validate_quantification(directory: Path, expected_version: str = '2.3.4',
                            expected_read_pairs: int | None = None) -> dict:
    names, tpm_total, reads_total = set(), 0.0, 0.0
    with (directory / 'quant.sf').open() as handle:
        table = csv.DictReader(handle, delimiter='\t')
        if table.fieldnames != ['Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads']:
            raise ValueError(f'Unexpected Salmon quant.sf schema: {directory}')
        for row in table:
            if None in row or not row['Name'] or row['Name'] in names:
                raise ValueError('Missing/duplicate transcript or malformed Salmon row')
            names.add(row['Name'])
            values = [float(row[key]) for key in ['Length', 'EffectiveLength', 'TPM', 'NumReads']]
            if any(not math.isfinite(value) or value < 0 for value in values) or values[0] <= 0:
                raise ValueError('Invalid Salmon numeric value')
            tpm_total += values[2]
            reads_total += values[3]
    with (directory / 'aux_info/meta_info.json').open() as handle:
        meta = json.load(handle)
    if not names or reads_total <= 0 or not math.isclose(tpm_total, 1e6, rel_tol=1e-4):
        raise ValueError('Empty/incomplete Salmon abundance table')
    if (str(meta.get('salmon_version')) != expected_version or
            type(meta.get('num_processed')) is not int or meta['num_processed'] <= 0):
        raise ValueError('Wrong Salmon version or no processed reads')
    if expected_read_pairs is not None:
        if type(expected_read_pairs) is not int or expected_read_pairs <= 0:
            raise ValueError('Expected paired-read count must be a positive integer')
        if meta.get('num_processed') != expected_read_pairs:
            raise ValueError('Salmon did not process every validated input read pair')
        if meta.get('em_converged') is not True or meta.get('quant_errors') != []:
            raise ValueError('Failed or unconverged Salmon abundance estimation')
    return {'transcripts': len(names), 'quant_sha256': digest(directory / 'quant.sf'),
            'meta_sha256': digest(directory / 'aux_info/meta_info.json')}


def require_provenance(path: Path, expected: dict) -> dict:
    if not path.is_file():
        raise ValueError(f'Unverified legacy cache: {path}. Move it aside and rebuild; existing files are preserved.')
    actual = json.loads(path.read_text())
    if actual.get('inputs') != expected:
        raise ValueError(f'Cache input/version mismatch: {path}. Move it aside and rebuild.')
    return actual


def write_provenance(path: Path, inputs: dict, outputs: dict) -> None:
    # Called only after successful completion and validation of all output files.
    import os
    import tempfile
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as handle:
            temporary = Path(handle.name)
            json.dump({'inputs': inputs, 'outputs': outputs}, handle, indent=2, sort_keys=True)
            handle.write('\n')
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def directory_fingerprint(directory: Path) -> dict:
    files = {str(path.relative_to(directory)): digest(path)
             for path in sorted(directory.rglob('*')) if path.is_file()}
    if not files:
        raise ValueError(f'Empty input directory: {directory}')
    return files
