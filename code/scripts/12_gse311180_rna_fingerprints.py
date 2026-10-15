#!/usr/bin/env python3
"""Guard canonical decompressed references and the full STAR index before RNA analysis.

The scientific producer is preserved byte for byte. Existing pre-guard results
require an explicit full-reference binding after their reconstruction is verified.
New runs bind references before producing any fingerprint; later changes fail.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import runpy
import sys
from input_integrity import digest
from rna_reference_integrity import bind_reference, verify_reference


def main():
    parser = argparse.ArgumentParser(description=__doc__, add_help=False)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bind-existing-verified-reference', action='store_true')
    parser.add_argument('--check-reference-only', action='store_true')
    if any(arg in ('-h', '--help') for arg in sys.argv[1:]):
        parser.print_help()
        print('Producer options: --threads N and --sample SAMPLE. Omit --sample for all 91 comparisons.')
        return
    args, forwarded = parser.parse_known_args()
    producer = Path(__file__).with_name('lib_rna_fingerprint_producer.py')
    references = verify_reference(args.reference.resolve())
    binding = bind_reference(args.output.resolve(), references,
                             allow_existing=args.bind_existing_verified_reference,
                             producer_sha256=digest(producer))
    print(json.dumps({'RNA_reference_artifact_checks': 'PASS',
                      'index_files_verified': len(references['index_files']),
                      'existing_run_bound_after_reconstruction': binding['existing_run_bound_after_reconstruction']}, indent=2), flush=True)
    if args.check_reference_only:
        return
    sys.argv = [str(producer), '--reference', str(args.reference), '--output', str(args.output), *forwarded]
    runpy.run_path(str(producer), run_name='__main__')


if __name__ == '__main__':
    main()
