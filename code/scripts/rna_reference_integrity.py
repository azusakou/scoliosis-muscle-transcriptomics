"""Full canonical RNA reference and immutable index binding for the RNA entry point."""
from __future__ import annotations
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import os
from input_integrity import digest, directory_fingerprint

PINNED = {
    'GRCh38.primary_assembly.genome.fa.gz': '3023c0705e83e86309271b09e25bdf678945fac98ce2fe5bb6b6e068f2225452',
    'gencode.v49.annotation.gtf.gz': 'd6e6fe0515c95b2a8cd36a853c1989cee9115c736c60237c56ae92b9daaaf7c4',
    'sites.hg38.rna.vcf.gz': '901cb1c396460a94a69b21f7f20e51fd2d74ed33ee0857343784bc123f4aa0fa',
}
PARAMETERS = {'versionGenome': '2.7.4a', 'genomeType': 'Full', 'genomeSAindexNbases': '14',
              'genomeChrBinNbits': '18', 'genomeSAsparseD': '1', 'genomeTransformType': 'None',
              'genomeTransformVCF': '-', 'sjdbOverhang': '149', 'sjdbFileChrStartEnd': '-',
              'sjdbGTFchrPrefix': '-', 'sjdbGTFfeatureExon': 'exon',
              'sjdbGTFtagExonParentTranscript': 'transcript_id', 'sjdbGTFtagExonParentGene': 'gene_id'}


def verify_reference(reference):
    compressed, uncompressed = {}, {}
    for name, expected in PINNED.items():
        path = reference / name
        actual = digest(path)
        if actual != expected:
            raise ValueError('Canonical reference checksum mismatch: ' + name)
        compressed[name] = actual
        if name.endswith(('.fa.gz', '.gtf.gz')):
            hashed = hashlib.sha256()
            with gzip.open(path, 'rb') as handle:
                for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
                    hashed.update(chunk)
            target = reference / name[:-3]
            target_hash = digest(target)
            if target_hash != hashed.hexdigest():
                raise ValueError('Uncompressed reference does not match canonical gzip: ' + target.name)
            uncompressed[target.name] = target_hash
    index = reference / 'STAR_2.7.11b'
    for name in ('Genome', 'SA', 'SAindex', 'genomeParameters.txt', 'chrNameLength.txt'):
        if not (index / name).is_file() or (index / name).stat().st_size == 0:
            raise ValueError('Incomplete STAR index: ' + name)
    parameters = {}
    for line in (index / 'genomeParameters.txt').read_text().splitlines():
        if line.startswith('#') or not line.strip():
            continue
        key, value = line.split(None, 1)
        parameters[key] = value.strip()
    for key, expected in PARAMETERS.items():
        if parameters.get(key) != expected:
            raise ValueError('Unexpected STAR index parameter: ' + key)
    for key, filename in [('genomeFastaFiles', 'GRCh38.primary_assembly.genome.fa'),
                          ('sjdbGTFfile', 'gencode.v49.annotation.gtf')]:
        value = parameters.get(key, '').split()
        if len(value) != 1 or Path(value[0]).name != filename:
            raise ValueError('Unexpected STAR indexed source: ' + key)
    return {'compressed': compressed, 'uncompressed': uncompressed,
            'index_files': directory_fingerprint(index), 'index_parameters': PARAMETERS}


def bind_reference(output, verified, *, allow_existing=False, producer_sha256):
    output.mkdir(parents=True, exist_ok=True)
    path = output / 'reference_artifact_binding.json'
    binding = {'references': verified, 'scientific_producer_sha256': producer_sha256}
    if path.exists():
        previous = json.loads(path.read_text())
        if previous['binding'] != binding:
            raise ValueError('RNA reference/index/producer changed; preserve existing results and use a new output directory')
        return previous
    existing = list(output.glob('*.provenance.json'))
    if existing and not allow_existing:
        raise ValueError('Older completed fingerprints lack full-reference binding; explicitly verify and bind that existing run, or use a new output directory')
    for stamp in existing:
        record = json.loads(stamp.read_text())
        if record['inputs']['recipe']['script_sha256'] != producer_sha256:
            raise ValueError('Existing fingerprint uses a different scientific producer')
    result = {'binding': binding, 'scope': 'Canonical source equivalence and full index/artifact identity; does not authenticate anatomical side',
              'existing_run_bound_after_reconstruction': bool(existing)}
    with tempfile.NamedTemporaryFile(mode='w', dir=output, encoding='utf-8', delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(result, handle, indent=2); handle.write('\n')
    os.replace(temporary, path)
    return result
