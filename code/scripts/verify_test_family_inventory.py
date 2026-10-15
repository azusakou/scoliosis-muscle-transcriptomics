#!/usr/bin/env python3
"""Reject omitted P-value artifacts in the retrospective analysis inventory.

An inventory block is not necessarily an independent testing family: archived
recipes, displayed subsets and repeated endpoint summaries overlap.
"""
import csv
import gzip
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def p_columns(columns):
    def raw_probability(column):
        key=column.lower().replace('.', '_')
        return (not key.startswith(('n_','reported_','p_reproduction_')) and
                not re.search(r'fdr|adj|min_attainable|minimum',key) and
                bool(re.search(r'(^|_)(pvalue|pval|p_value|p_val|p)(_|$)|_tail$',key)))
    return [c for c in columns if raw_probability(c)]


def validate(root, records):
    names = [r['family'] for r in records]
    if len(set(names)) != len(names):
        raise ValueError('Duplicate inventory block')
    listed = {r['source'] for r in records}
    for r in records:
        if not (root / r['source']).is_file() or not r['role'] or not r['correction']:
            raise ValueError('Missing inventory source/interpretation')
    missing = []
    artifacts = 0
    for directory in ['results/tables', 'analysis', 'results/figures/source_data', 'analysis/tables', 'analysis/bulk_gate_inputs']:
        parent = root / directory
        paths = parent.rglob('*.tsv*') if directory == 'results/tables' else parent.glob('*.tsv*')
        for path in paths:
            with (gzip.open(path, 'rt') if path.name.endswith('.gz') else path.open()) as f:
                reader = csv.reader(f, delimiter='\t')
                columns = next(reader)
            if p_columns(columns):
                artifacts += 1
                relative = path.relative_to(root).as_posix()
                if relative not in listed:
                    missing.append(relative)
    if missing:
        raise ValueError('Omitted P-value artifacts: ' + '; '.join(sorted(missing)))
    return dict(status='PASS',inventory_blocks=len(records),P_value_artifacts_covered=artifacts,
                scope='Source-artifact coverage, including duplicate summaries; no independence or study-wide FDR claimed')


if __name__ == '__main__':
    with (ROOT / 'analysis/bulk_test_family_inventory.tsv').open() as f:
        records = list(csv.DictReader(f, delimiter='\t'))
    print(json.dumps(validate(ROOT, records), indent=2))
