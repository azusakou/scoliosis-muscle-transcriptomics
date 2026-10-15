# Supplementary RNA profile pairing check

This is separate from Salmon gene quantification. It compares all 14 libraries
(91 pairs) against the seven metadata-defined pairs. The interpretation gates
were recorded before extraction in
`metadata/audit/rna_identity_consistency_criteria.json`. Agreement supports
RNA profile pairing; it does not identify a person, prove anatomical side,
resolve the deposited processed-matrix anomaly, or supply missing clinical data.
Low coverage is inconclusive. Do not tune gates after observing the results.

Create the verified `environment.identity.linux-64.explicit.txt` lock in its own prefix (`environment_identity.yml` is the bootstrap alternative) and place its `bin` directory
on PATH. This 30-artifact tools lock does not include Python. Set
`AIS_ANALYSIS_PYTHON` to the Python3.12 interpreter in the installed analysis
environment; do not rely on a system Python selected after changing PATH.
Stage 08 must have completed full paired FASTQ validation and written
the source-bound per-sample provenance files. Keep the same `AIS_RAW_ROOT`.
The reference index requires substantial RAM/disk; the actual run used 96 GB
as the STAR build limit, 18 build threads and separate SSD scratch. Each human
alignment needs about 32 GB for its genome index. The public script aligns one
sample at a time by default; the final uniform run uses three concurrent eight-thread jobs.
The seed window cap is 20; STAR writes unsorted BAM and samtools1.24 sorts
with 512 MiB per sorting thread. Reserve about 32 GB per simultaneous STAR
job, plus sorting/I/O buffers and operating-system space.

Prepare the public references below in an isolated directory. These commands
are for Bash on Linux; use full paths in the two task-specific variables.
Downloads must match the pinned SHA-256 values before decompression or indexing.
The reference directory must be new or contain only independently verified
matching inputs. Preserve incomplete/unverified index directories separately.

```bash
set -euo pipefail
export AIS_RNA_REFERENCE=/absolute/path/to/new/rna-reference
export AIS_RNA_OUTPUT=/absolute/path/to/new/rna-output
mkdir -p "$AIS_RNA_REFERENCE"
curl --fail --location --retry 8 \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz \
  -o "$AIS_RNA_REFERENCE/GRCh38.primary_assembly.genome.fa.gz"
curl --fail --location --retry 8 \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz \
  -o "$AIS_RNA_REFERENCE/gencode.v49.annotation.gtf.gz"
curl --fail --location --retry 8 \
  https://github.com/brentp/somalier/files/4566475/sites.hg38.rna.vcf.gz \
  -o "$AIS_RNA_REFERENCE/sites.hg38.rna.vcf.gz"
"$AIS_ANALYSIS_PYTHON" - <<'PY'
import gzip, hashlib, json, os
from pathlib import Path
r = Path(os.environ['AIS_RNA_REFERENCE'])
sources = [
 ('GRCh38.primary_assembly.genome.fa.gz',
  'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz',
  '3023c0705e83e86309271b09e25bdf678945fac98ce2fe5bb6b6e068f2225452'),
 ('gencode.v49.annotation.gtf.gz',
  'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz',
  'd6e6fe0515c95b2a8cd36a853c1989cee9115c736c60237c56ae92b9daaaf7c4'),
 ('sites.hg38.rna.vcf.gz',
  'https://github.com/brentp/somalier/files/4566475/sites.hg38.rna.vcf.gz',
  '901cb1c396460a94a69b21f7f20e51fd2d74ed33ee0857343784bc123f4aa0fa')]
manifest = []
for name, url, expected in sources:
    p = r / name
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''): h.update(b)
    if h.hexdigest() != expected: raise ValueError('Reference SHA mismatch: ' + name)
    with gzip.open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''): pass
    manifest.append({'file': name, 'url': url, 'bytes': p.stat().st_size, 'sha256': expected})
(r / 'reference_sha256_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
PY
gzip -dc "$AIS_RNA_REFERENCE/GRCh38.primary_assembly.genome.fa.gz" \
  > "$AIS_RNA_REFERENCE/GRCh38.primary_assembly.genome.fa.tmp"
mv "$AIS_RNA_REFERENCE/GRCh38.primary_assembly.genome.fa.tmp" \
  "$AIS_RNA_REFERENCE/GRCh38.primary_assembly.genome.fa"
gzip -dc "$AIS_RNA_REFERENCE/gencode.v49.annotation.gtf.gz" \
  > "$AIS_RNA_REFERENCE/gencode.v49.annotation.gtf.tmp"
mv "$AIS_RNA_REFERENCE/gencode.v49.annotation.gtf.tmp" \
  "$AIS_RNA_REFERENCE/gencode.v49.annotation.gtf"
if [[ -e "$AIS_RNA_REFERENCE/STAR_2.7.11b" ]]; then
  printf 'Preserve and move the existing index before rebuilding.\n' >&2
  exit 1
fi
[[ "$(STAR --version)" == "2.7.11b" ]]
mkdir -p "$AIS_RNA_REFERENCE/STAR_2.7.11b"
STAR --runMode genomeGenerate --runThreadN 18 \
  --limitGenomeGenerateRAM 96000000000 \
  --genomeDir "$AIS_RNA_REFERENCE/STAR_2.7.11b" \
  --genomeFastaFiles "$AIS_RNA_REFERENCE/GRCh38.primary_assembly.genome.fa" \
  --sjdbGTFfile "$AIS_RNA_REFERENCE/gencode.v49.annotation.gtf" \
  --sjdbOverhang 149 \
  --outFileNamePrefix "$AIS_RNA_REFERENCE/STAR_build_"
"$AIS_ANALYSIS_PYTHON" scripts/12_gse311180_rna_fingerprints.py \
  --reference "$AIS_RNA_REFERENCE" --output "$AIS_RNA_OUTPUT" --threads 8
"$AIS_ANALYSIS_PYTHON" scripts/13_export_rna_identity_audit.py \
  --identity-output "$AIS_RNA_OUTPUT" \
  --raw-root "$AIS_RAW_ROOT" \
  --report "$AIS_RNA_OUTPUT/completed_raw_pairing_export.json"
"$AIS_ANALYSIS_PYTHON" scripts/verify_rna_identity.py \
  --pairs analysis/audit_gse311180_rna_all_pairs.tsv \
  --manifest metadata/audit/gse311180_requant_manifest.tsv \
  --criteria metadata/audit/rna_identity_consistency_criteria.json \
  --expected-pairs analysis/audit_gse311180_rna_expected_pairs.tsv --mapping-qc analysis/audit_gse311180_raw_RNA_mapping_QC.tsv \
  --output "$AIS_RNA_OUTPUT/independent_completed_pair_assessment.json"
Rscript analysis/scripts/figure_rna_pairing.R

"$AIS_ANALYSIS_PYTHON" scripts/verify_rna_identity.py \
  --pairs "$AIS_RNA_OUTPUT/rna_identity.pairs.tsv" \
  --manifest metadata/audit/gse311180_requant_manifest.tsv \
  --criteria metadata/audit/rna_identity_consistency_criteria.json \
  --output "$AIS_RNA_OUTPUT/independent_pair_assessment.json"
```

The all-sample script verifies existing FASTQ/BAM/fingerprint checksums before
reuse, preserves successful per-sample checkpoints, and rejects unverified
partial caches. `--sample` creates a one-sample checkpoint; omit it to compare
all samples. Alignment is two-pass STAR with `--seedPerWindowNmax 20` and unsorted BAM
output, followed by `samtools sort -@ 8 -m 512M` and coordinate-BAM integrity
checks. These fixed settings apply to every sample. The initial default-cap
attempt encountered recursive `stitchWindowAligns` slowdown; the
[STAR maintainer documents this issue](https://github.com/alexdobin/STAR/issues/1414#issuecomment-991980437).
The lower cap resolved the affected full-read mapping trial, but can change
mapping sensitivity. The trial then required more than 12.23 GB for STAR
sorting. To bound memory across libraries without changing mapped alignments,
sorting was delegated to samtools after STAR exits. Neither adjustment changes
the Salmon expression analysis or the recorded pairing interpretation gates.
Old attempts are preserved separately and are excluded from the final uniform
recipe. Failed stages are preserved for diagnosis and cannot be silently reused. Somalier uses
the pinned hg38 RNA-specific sites, minimum depth 7 and heterozygote allele
balance 0.2. Read the complete pair matrix and individual mapping/coverage
metrics before interpretation. The independent validator checks all 91 pairs,
relatedness arithmetic, callable-site coverage, IBS0 fraction and reciprocal
unique best matches. A structural PASS is distinct from a consistency verdict.
The validator recomputes relatedness from the integer genotype counts before
applying the gates and nearest-match comparisons. Somalier's three-decimal
display value is checked against that formula but is not used for thresholding.

The portable exporter checks all-library uniform recipes, completed Salmon/source
provenance, full-read counts, native fingerprint hashes and mapping/coverage QC
before exporting four aggregate tables. It exports no genotype sketches or
variant calls. Run stage08 and the all-sample stage12 source validation first.

Native tools remain under their own licences. No BAMs, read-level SNP calls,
human genotype sketches, software binaries or genome references are redistributed
in the code release.


The public stage12 entry point now verifies that both decompressed genome and
GTF match their pinned gzip sources, checks the fixed STAR index parameters,
and hashes every index file. New output directories bind these full reference
artifacts before processing. Changed references, index bytes or scientific
producer subsequently stop cache reuse. The scientific producer is preserved
byte for byte in `lib_rna_fingerprint_producer.py`; the guarded entry point
executes that same producer.

An older already verified reconstruction needs the explicit
`--bind-existing-verified-reference` flag on its first guarded invocation.
This creates a binding after that existing reconstruction and records that
fact; it does not retroactively claim a prior full-index hash record. The
current study's canonical decompressed references/index were checked, and its
completed DR50 cache was tested through this guarded entry point. The new
regression tests reject stale decompressed sources, changed index bytes,
changed index parameters and unapproved/different-producer legacy bindings.

The supplied metadata reference-binding record describes this particular actual
index build. STAR parameter metadata contain build paths, so rebuilt indexes
need not be byte-identical across locations. Bind a newly built index in a new
output directory; do not copy the supplied study binding into another build.
The compressed/decompressed canonical sequence and annotation bytes remain
fixed, and index parameter checks apply independently of build path.

## Post hoc investigation of STB-084 and STB-100

Keep the primary cap20 directory and recorded gates unchanged. All settings below
are diagnostic sensitivities; a higher-depth score does not overturn a primary
failure or authenticate DNA identity/anatomical side. Use a new private output
root and sufficient RAM (about 32 GiB per simultaneous STAR job plus buffers).
These commands start from the actual completed 14-library primary output.

```bash
export AIS_RNA_DIAGNOSTICS=/absolute/path/to/new/private/rna-diagnostics
mkdir -p "$AIS_RNA_DIAGNOSTICS"
"$AIS_ANALYSIS_PYTHON" scripts/14_rna_pairing_failure_diagnostics.py \
  --identity "$AIS_RNA_OUTPUT" --sites "$AIS_RNA_REFERENCE/sites.hg38.rna.vcf.gz" \
  --default-pairs "$AIS_RNA_OUTPUT/rna_identity.pairs.tsv" \
  --output "$AIS_RNA_DIAGNOSTICS/primary_counts"
"$AIS_ANALYSIS_PYTHON" scripts/20_rna_native_diagnostic_grid.py \
  --identity "$AIS_RNA_OUTPUT" --sites "$AIS_RNA_REFERENCE/sites.hg38.rna.vcf.gz" \
  --output "$AIS_RNA_DIAGNOSTICS/native_grid"
"$AIS_ANALYSIS_PYTHON" scripts/17_rna_dropout_context.py \
  --identity "$AIS_RNA_OUTPUT" --sites "$AIS_RNA_REFERENCE/sites.hg38.rna.vcf.gz" \
  --output "$AIS_RNA_DIAGNOSTICS/primary_counts"
"$AIS_ANALYSIS_PYTHON" scripts/16_rna_unique_read_pileup_diagnostics.py \
  --identity "$AIS_RNA_OUTPUT" --reference "$AIS_RNA_REFERENCE" \
  --output "$AIS_RNA_DIAGNOSTICS/strict_reads" --paired-overlap remove --jobs 2 --view-threads 2
"$AIS_ANALYSIS_PYTHON" scripts/19_rna_read_filter_transitions.py \
  --primary "$AIS_RNA_OUTPUT" --strict "$AIS_RNA_DIAGNOSTICS/strict_reads" \
  --output "$AIS_RNA_DIAGNOSTICS/strict_reads"
"$AIS_ANALYSIS_PYTHON" scripts/23_rna_fragment_coordinate_diversity.py \
  --identity "$AIS_RNA_OUTPUT" --sites "$AIS_RNA_REFERENCE/sites.hg38.rna.vcf.gz" \
  --output "$AIS_RNA_DIAGNOSTICS/fragment_coordinates"
"$AIS_ANALYSIS_PYTHON" scripts/15_rna_pairing_cap50_sensitivity.py \
  --raw "$AIS_RAW_ROOT" --reference "$AIS_RNA_REFERENCE" \
  --identity "$AIS_RNA_OUTPUT" --output "$AIS_RNA_DIAGNOSTICS/cap50" \
  --jobs 2 --threads 8 --timeout 10800
"$AIS_ANALYSIS_PYTHON" scripts/18_export_rna_cap50_diagnostics.py \
  --primary "$AIS_RNA_OUTPUT" --completed "$AIS_RNA_DIAGNOSTICS/cap50" --output "$AIS_RNA_DIAGNOSTICS/cap50"
"$AIS_ANALYSIS_PYTHON" scripts/21_export_rna_cap50_attempt_status.py \
  --runs "$AIS_RNA_DIAGNOSTICS/cap50" --primary "$AIS_RNA_OUTPUT" \
  --output "$AIS_RNA_DIAGNOSTICS/cap50"
# Copy only verified audit_*.tsv aggregates from the private subdirectories
# into analysis; never copy BAMs/sketches/private marker-count files.
"$AIS_ANALYSIS_PYTHON" scripts/verify_rna_failure_diagnostics.py
Rscript analysis/scripts/figure_rna_failure_diagnostics.R
```

Script 15's default sample
subset contains both nonpassing pairs and the passing STB-085 control; arbitrary
subsets must contain complete metadata pairs. The source-only guard mode
`--validate-inputs-only` checks whole FASTQ bytes, pinned references/full index,
primary provenance and tool versions without performing a new alignment. Such a
mode is not evidence that a sensitivity alignment completed.

The format2 decoder uses lexically sorted autosomes and alphabetically ordered
allele counts, rather than assuming reference/alternate count order. The native
v0.3.5 third-base ratio is rounded to float32 before its 0.1 cutoff; adversarial
boundary fixtures and actual 1,092 native comparisons check this detail. All
primary integer counts reproduce. The strict pileup uses MAPQ60/baseQ20 and uncapped depth, with BAQ disabled.
The default `--paired-overlap remove` uses samtools paired-overlap quality
adjustment; `--paired-overlap retain` adds `-x` and counts both overlapping mate
bases. Both modes are reported separately, and neither is presented as exactly
the native Somalier/Hileup counting implementation. Original
BAMs are not PCR duplicate-marked; filtering duplicate flags is not deduplication.
The native CHARR/paired-likelihood tables are retained only as uncalibrated RNA
proxies and must not be interpreted as contamination percentages. See
`analysis/RNA_PAIRING_FAILURE_DIAGNOSTIC_REPORT.md` for the observed loss,
coverage, same-marker transitions, complete grid and alignment sensitivity.

The earlier cap50 alignment producer is preserved byte for byte as a non-entry-point
execution record in `scripts/execution_records/rna_cap50_alignment_producer_20261013.txt`.
Use the guarded stage 15 command, not that archived record. Its STAR/sort/extract
scientific commands are retained; the guarded entry additionally validates full
canonical references/indexes, provenance/tool versions and complete sample-pair
arguments. Execution records/receipt hashes distinguish the two producer versions
without asserting that later guard checks were recorded before earlier alignments.
The separately executed guarded producer is also preserved in
`scripts/execution_records/rna_cap50_guarded_alignment_producer_20261013.txt`.
The current entry point additionally reports the actual selected-library count
in its completion label; archived receipts remain the authoritative record of
which libraries completed.

An initial diagnostic proof label incorrectly called the overlap-retaining
pileup the original/native counting convention. Its actual command used `-x`;
the retained numbers are preserved under explicitly named overlap-retained
source tables. The native source delegates to Hileup, which handles overlapping
mates. The additional overlap-adjusted run addresses this distinction. The
archived execution text and correction preserve the audit trail rather than
rewriting the historical producer hash.

For the separate overlap-retaining check, run stage16 in another new private
output directory with `--paired-overlap retain`, then stage19 against that
same directory with `--overlap-retained`. Copy only its explicitly suffixed
aggregate TSVs. The default adjustment and retention modes use the same
MAPQ/baseQ/site filters; neither is asserted to duplicate native extraction.


The current cap50 runner starts native stages in an isolated POSIX process
group and stops the entire group on timeout/failure. This prevents a STAR shell
wrapper from leaving a native child writing into a failed stage. The regression
check observes that a spawned child's file-writing heartbeat stops on timeout.
Do not treat a leftover BAM or continued legacy child as a verified completion;
require the complete source-bound receipt and full input/output checks.

Stage23 checks all original discordant loci in all14 BAMs. It parses spliced
CIGAR/base-quality records, requires proper pairs and excludes supplementary
alignments, collapses overlapping mates by QNAME/template coordinates and
reports only directional aggregates. Coordinate diversity is not UMI-based
molecule counting, PCR duplicate measurement or DNA genotyping. Private QNAMEs,
marker alleles and coordinates must remain on the compute host.

Two legacy cap50 wrappers timed out while their native STAR children continued.
Their FAILED_PRESERVED receipts remain unchanged. After full native completion,
stage22 verified raw/reference/index bytes, exact STAR BAM header command/RG,
ALL DONE in Log.out, successful stdout, full input counts and a complete BAM,
then sorted/indexed/extracted separately with the original BAM hash unchanged.
The initial recovery guard checked the wrong log and failed before output;
that attempt is retained privately. This is recovery of a COMPLETE alignment,
never acceptance of a partial alignment. Use stage22 only with its explicitly
supported DR118/DR148 legacy source-bound receipt and a new namespace. Current
stage15 kills the entire process group on timeout, so it needs no such recovery.

When complete pair outputs come from separate verified namespaces, create a
private master directory with sample-name links to exactly one PASS directory
per sample. Generate the native comparison afresh in a new output directory:

```bash
"$AIS_ANALYSIS_PYTHON" scripts/18_export_rna_cap50_diagnostics.py \
  --primary "$AIS_RNA_OUTPUT" --completed /private/complete/master \
  --output "$AIS_RNA_DIAGNOSTICS/cap50_export" --generate-native-pairs \
  --sites "$AIS_RNA_REFERENCE/sites.hg38.rna.vcf.gz"
"$AIS_ANALYSIS_PYTHON" scripts/21_export_rna_cap50_attempt_status.py \
  --runs /private/complete/master --primary "$AIS_RNA_OUTPUT" \
  --output "$AIS_RNA_DIAGNOSTICS/cap50_attempts"
```

Stage18 requires full BAM/fingerprint hashes and independently reproduces every
native integer comparison. Stage21 collects source-bound terminal metadata,
including explicit legacy recovery status; it does not replace BAM validation.
Do not pass both a legacy failure and its new recovery as competing receipts.
Executed sources and sanitized scope/table hashes are recorded in
metadata/audit/rna_failure_diagnostics_execution.json.
