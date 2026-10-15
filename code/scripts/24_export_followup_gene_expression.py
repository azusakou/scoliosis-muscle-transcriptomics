#!/usr/bin/env python3
"""Export all mapped gene TPM/read counts from the verified final Salmon run.

This is an expression audit, never a genotype/side identity assay. Older
quantification archives are rejected against the released source hashes.
"""
import argparse
import csv
import gzip
import io
import json
from decimal import Decimal, InvalidOperation
from pathlib import Path
from input_integrity import digest, validate_quantification

ROOT = Path(__file__).resolve().parents[1]


def aggregate_gene_rows(rows, mapping):
    # Quantification files print decimal values. Exact decimal summation avoids
    # creating unequal binary last bits for mathematically tied gene TPMs.
    genes = {}; mapped = total = 0
    for row in rows:
        total += 1; gene = mapping.get(row['Name'])
        if not gene:
            continue
        try:
            values = [Decimal(row[k]) for k in ('TPM', 'NumReads')]
        except InvalidOperation as error:
            raise ValueError('Malformed quantification number') from error
        if any(not v.is_finite() or v<0 for v in values):
            raise ValueError('Invalid quantification number')
        mapped += 1; value = genes.setdefault(gene, [Decimal(0), Decimal(0)])
        for i in (0, 1):
            value[i] += values[i]
    return genes, mapped, total


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--quant-root', type=Path, required=True)
    p.add_argument('--tx2gene', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    source = ROOT / 'analysis/audit_gse311180_completed_raw_source_sha256.tsv'
    qc = ROOT / 'analysis/audit_gse311180_raw_RNA_mapping_QC.tsv'
    reference = ROOT / 'metadata/audit/gencode_v49_reference.tsv'
    with source.open() as f:
        records = list(csv.DictReader(f, delimiter='\t'))
    with qc.open() as f:
        counts = {r['sample_id']: int(r['validated_read_pairs']) for r in csv.DictReader(f, delimiter='\t')}
    with reference.open() as f:
        ref = {r['field']: r['value'] for r in csv.DictReader(f, delimiter='\t')}
    if digest(a.tx2gene) != ref['tx2gene_sha256']:
        raise ValueError('Transcript-to-gene reference differs from the final source')
    with a.tx2gene.open() as f:
        mapping = {r['transcript_id']: r['gene_name'] for r in csv.DictReader(f, delimiter='\t')}
    vectors = {}; proof = []
    for record in records:
        name = record['sample_id']; directory = a.quant_root / name
        validation = validate_quantification(directory, expected_read_pairs=counts[name])
        if (validation['quant_sha256'] != record['quant_sha256'] or
                validation['meta_sha256'] != record['quant_metadata_sha256']):
            raise ValueError('Not the final verified quantification: ' + name)
        with (directory / 'quant.sf').open() as f:
            genes, mapped, total = aggregate_gene_rows(csv.DictReader(f, delimiter='\t'), mapping)
        if mapped / total < .95:
            raise ValueError('Insufficient transcript-to-gene mapping')
        vectors[name] = genes
        proof.append(dict(sample_id=name, mapped_transcripts=mapped,
                          total_transcripts=total, **validation))
    if len(vectors) != 14:
        raise ValueError('Require all 14 libraries')
    a.output.parent.mkdir(parents=True, exist_ok=True)
    names = sorted(vectors); genes = sorted(set().union(*(set(v) for v in vectors.values())))
    with a.output.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding='utf-8', newline='') as text:
                w = csv.writer(text, delimiter='\t', lineterminator='\n')
                w.writerow(['gene'] + [n+'_'+kind for n in names for kind in ['TPM', 'NumReads']])
                for gene in genes:
                    w.writerow([gene] + [value for n in names for value in vectors[n].get(gene, [0., 0.])])
    a.output.with_suffix(a.output.suffix+'.proof.json').write_text(json.dumps(dict(
        status='PASS', source_sha256=digest(source), tx2gene_sha256=digest(a.tx2gene),
        producer_sha256=digest(Path(__file__)), samples=14, gene_rows=len(genes),
        output_sha256=digest(a.output), quantification_checks=proof,
        scope='Unfiltered gene-summed expression only; no side/DNA authentication'), indent=2)+'\n')
    print('FINAL_SOURCE_GENE_EXPRESSION_EXPORT_PASS', len(genes), flush=True)


if __name__ == '__main__':
    main()
