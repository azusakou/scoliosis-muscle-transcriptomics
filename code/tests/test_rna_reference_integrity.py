"""Regression tests for previously unbound decompressed references/index caches."""
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from rna_reference_integrity import PINNED, PARAMETERS, bind_reference, verify_reference


class ReferenceGuards(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        pins = {}
        for name in PINNED:
            payload = b'>toy\nACGT\n' if '.fa.' in name else b'toy source\n'
            path = self.root / name
            path.write_bytes(gzip.compress(payload, mtime=0))
            pins[name] = hashlib.sha256(path.read_bytes()).hexdigest()
            if name.endswith(('.fa.gz', '.gtf.gz')):
                (self.root / name[:-3]).write_bytes(payload)
        self.index = self.root / 'STAR_2.7.11b'
        self.index.mkdir()
        for name in ('Genome', 'SA', 'SAindex', 'chrNameLength.txt'):
            (self.index / name).write_text('toy index\n')
        params = dict(PARAMETERS, genomeFastaFiles=str(self.root / 'GRCh38.primary_assembly.genome.fa'),
                      sjdbGTFfile=str(self.root / 'gencode.v49.annotation.gtf'))
        (self.index / 'genomeParameters.txt').write_text(''.join(k+'\t'+v+'\n' for k,v in params.items()))
        mocked = patch('rna_reference_integrity.PINNED', pins)
        mocked.start(); self.addCleanup(mocked.stop)
        self.output = self.root / 'outputs'
        self.producer = 'a' * 64

    def test_valid_and_unchanged_binding(self):
        verified = verify_reference(self.root)
        first = bind_reference(self.output, verified, producer_sha256=self.producer)
        self.assertFalse(first['existing_run_bound_after_reconstruction'])
        self.assertEqual(first, bind_reference(self.output, verified, producer_sha256=self.producer))

    def test_correct_gzip_cannot_hide_stale_uncompressed_genome(self):
        (self.root / 'GRCh38.primary_assembly.genome.fa').write_text('WRONG GENOME\n')
        with self.assertRaisesRegex(ValueError, 'Uncompressed reference'):
            verify_reference(self.root)

    def test_changed_index_rejects_existing_cache_binding(self):
        bind_reference(self.output, verify_reference(self.root), producer_sha256=self.producer)
        (self.index / 'SA').write_text('different index bytes\n')
        with self.assertRaisesRegex(ValueError, 'reference/index/producer changed'):
            bind_reference(self.output, verify_reference(self.root), producer_sha256=self.producer)

    def test_wrong_index_parameters_rejected(self):
        p = self.index / 'genomeParameters.txt'
        p.write_text(p.read_text().replace('sjdbOverhang\t149', 'sjdbOverhang\t148'))
        with self.assertRaisesRegex(ValueError, 'index parameter'):
            verify_reference(self.root)

    def test_legacy_binding_is_explicit_and_same_producer_only(self):
        self.output.mkdir()
        stamp = self.output / 'toy.provenance.json'
        stamp.write_text(json.dumps({'inputs': {'recipe': {'script_sha256': self.producer}}}))
        verified = verify_reference(self.root)
        with self.assertRaisesRegex(ValueError, 'lack full-reference binding'):
            bind_reference(self.output, verified, producer_sha256=self.producer)
        stamp.write_text(json.dumps({'inputs': {'recipe': {'script_sha256': 'b'*64}}}))
        with self.assertRaisesRegex(ValueError, 'different scientific producer'):
            bind_reference(self.output, verified, allow_existing=True, producer_sha256=self.producer)
        stamp.write_text(json.dumps({'inputs': {'recipe': {'script_sha256': self.producer}}}))
        result = bind_reference(self.output, verified, allow_existing=True, producer_sha256=self.producer)
        self.assertTrue(result['existing_run_bound_after_reconstruction'])


if __name__ == '__main__':
    unittest.main()
