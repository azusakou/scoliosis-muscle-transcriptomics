#!/usr/bin/env python3
"""Supplementary RNA fingerprint check; anatomical side remains metadata-defined.

Requires validated paired FASTQs from stage 08, STAR 2.7.11b index made with
GENCODE49 primary GRCh38 and Somalier 0.3.5 RNA-specific hg38 sites.
Runs one alignment at a time. Uses a fixed seed-window cap of 20 after
observed default-50 recursive stitching slowdown. Does not change the
Salmon transcriptome analysis; interpretation thresholds remain predeclared.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import subprocess
import tempfile
from pathlib import Path

from input_integrity import digest, require_provenance, write_provenance, validate_pairs
from project_paths import raw_data_root, executable

ROOT = Path(__file__).resolve().parents[1]
FASTA_SHA = "3023c0705e83e86309271b09e25bdf678945fac98ce2fe5bb6b6e068f2225452"
SITES_SHA = "901cb1c396460a94a69b21f7f20e51fd2d74ed33ee0857343784bc123f4aa0fa"


def run(argv, logfile):
    with logfile.open("w") as handle:
        subprocess.run(argv, stdout=handle, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--sample", help="Optional one-sample checkpoint; relate runs only with all 14 digests")
    args = parser.parse_args()
    assert args.threads > 0
    reference, output = args.reference.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    fastq = raw_data_root() / "bulk_gse311180/fastq"
    with (ROOT / "metadata/audit/gse311180_requant_manifest.tsv").open() as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    validate_pairs(manifest, expected_pairs=7)
    if args.sample and args.sample not in {r["sample_id"] for r in manifest}:
        raise ValueError("Unknown sample")
    star = executable("STAR", "AIS_STAR")
    samtools = executable("samtools", "AIS_SAMTOOLS")
    somalier = executable("somalier", "AIS_SOMALIER")
    samtools_version = subprocess.check_output([samtools, "--version"], text=True).splitlines()[0]
    assert samtools_version == "samtools 1.24", samtools_version
    star_version = subprocess.check_output([star, "--version"], text=True).strip()
    assert star_version == "2.7.11b", star_version
    help_result = subprocess.run([somalier, "--help"], capture_output=True, text=True)
    assert "somalier version: 0.3.5" in help_result.stdout + help_result.stderr
    genome = reference / "GRCh38.primary_assembly.genome.fa"
    sites = reference / "sites.hg38.rna.vcf.gz"
    index = reference / "STAR_2.7.11b"
    assert digest(genome.with_suffix(genome.suffix + ".gz")) == FASTA_SHA
    assert digest(sites) == SITES_SHA
    for filename in ("Genome", "SA", "SAindex", "genomeParameters.txt"):
        assert (index / filename).is_file(), f"Incomplete STAR index: {filename}"
    reference_manifest = reference / "reference_sha256_manifest.json"
    reference_hash = digest(reference_manifest)
    if not genome.with_suffix(genome.suffix + ".fai").exists():
        run([samtools, "faidx", str(genome)], output / "reference_faidx.log")
    recipe = {"STAR": star_version, "Somalier": "0.3.5", "reference_manifest_sha256": reference_hash,
              "script_sha256": digest(Path(__file__)),
              "sites_sha256": SITES_SHA, "twopassMode": "Basic", "unique_mapq": 60,
              "minimum_genotype_depth": 7, "minimum_heterozygote_allele_balance": .2,
              "seedPerWindowNmax": 20, "STAR_BAM_output": "Unsorted",
              "samtools": "1.24", "samtools_sort_memory_per_thread": "512M",
              "alignment_threads": args.threads,
              "interpretation": "Genotype-like RNA pairing check; cannot authenticate anatomical side or clinical identity"}
    (output / "method_configuration.json").write_text(json.dumps(recipe, indent=2) + "\n")
    fingerprints = output / "fingerprints"
    fingerprints.mkdir(exist_ok=True)
    for row in manifest:
        sample = row["sample_id"]
        if args.sample and sample != args.sample:
            continue
        fq_stamp = fastq / f"{sample}.provenance.json"
        fq = json.loads(fq_stamp.read_text())
        assert fq["inputs"]["sample_id"] == sample and fq["inputs"]["sra_run"] == row["sra_run"]
        mates = [fastq / f"{sample}_{mate}.fastq" for mate in (1, 2)]
        for mate, path in enumerate(mates, 1):
            assert digest(path) == fq["outputs"][f"mate{mate}_sha256"]
        inputs = {"sample": sample, "fastq_provenance_sha256": digest(fq_stamp), "recipe": recipe}
        stamp = output / f"{sample}.provenance.json"
        bam = output / sample / "Aligned.sortedByCoord.out.bam"
        fingerprint = fingerprints / f"{sample}.somalier"
        if stamp.exists():
            saved = require_provenance(stamp, inputs)
            assert saved["outputs"]["bam_sha256"] == digest(bam)
            assert saved["outputs"]["fingerprint_sha256"] == digest(fingerprint)
            subprocess.run([samtools, "quickcheck", "-v", str(bam)], check=True)
            print(sample, "VERIFIED_EXISTING", flush=True)
            continue
        if (output / sample).exists() or fingerprint.exists():
            raise ValueError(f"Unverified partial alignment/cache for {sample}; preserve it separately before retrying")
        stage = Path(tempfile.mkdtemp(prefix=f"{sample}_stage_", dir=output))
        try:
            run([star, "--runThreadN", str(args.threads), "--genomeDir", str(index),
                 "--readFilesIn", *(str(p) for p in mates), "--twopassMode", "Basic",
                 "--outSAMtype", "BAM", "Unsorted", "--outSAMmapqUnique", "60",
                 "--seedPerWindowNmax", "20",
                 "--outSAMattrRGline", f"ID:{sample}", f"SM:{sample}", "PL:ILLUMINA",
                 "--outFileNamePrefix", str(stage) + "/"], output / f"{sample}.alignment.log")
            staged_bam = stage / bam.name
            unsorted_bam = stage / "Aligned.out.bam"
            run([samtools, "sort", "-@", str(args.threads), "-m", "512M",
                 "-o", str(staged_bam), str(unsorted_bam)], output / f"{sample}.sort.log")
            subprocess.run([samtools, "quickcheck", "-v", str(staged_bam)], check=True)
            run([samtools, "index", "-@", str(args.threads), str(staged_bam)], output / f"{sample}.index.log")
            run([samtools, "flagstat", str(staged_bam)], stage / "samtools_flagstat.txt")
            run([somalier, "extract", "--sites", str(sites), "--fasta", str(genome),
                 "--out-dir", str(stage), str(staged_bam)], output / f"{sample}.extract.log")
            staged_fingerprint = stage / fingerprint.name
            assert staged_fingerprint.is_file() and staged_fingerprint.stat().st_size > 0
            outputs = {"bam_sha256": digest(staged_bam), "fingerprint_sha256": digest(staged_fingerprint),
                       "STAR_log_sha256": digest(stage / "Log.final.out")}
            staged_fingerprint.replace(fingerprint)
            # This generated intermediate is redundant after validated sorting.
            unsorted_bam.unlink()
            # Move completed files atomically on the same filesystem.
            retained = stage / "retained"
            retained.mkdir()
            for file in list(stage.iterdir()):
                if file.is_file():
                    file.replace(retained / file.name)
            retained.replace(output / sample)
        except BaseException:
            # Preserve logs and partial outputs for diagnosis; never reuse them.
            if stage.exists():
                stage.rename(stage.with_name(stage.name + ".failed"))
            raise
        else:
            import shutil
            shutil.rmtree(stage)
        write_provenance(stamp, inputs, outputs)
        print(sample, "ALIGNMENT_AND_FINGERPRINT_PASS", flush=True)
    if not args.sample and all((output / f'{r["sample_id"]}.provenance.json').is_file() for r in manifest):
        groups = output / "expected_pairs.txt"
        patients = sorted({r["patient_id"] for r in manifest})
        groups.write_text("".join(",".join(r["sample_id"] for r in manifest if r["patient_id"] == patient) + "\n" for patient in patients))
        env = dict(os.environ, SOMALIER_REPORT_ALL_PAIRS="1")
        with (output / "relate.log").open("w") as handle:
            subprocess.run([somalier, "relate", "--sites", str(sites), "--groups", str(groups),
                            "--min-depth", "7", "--min-ab", "0.2", "--output-prefix", str(output / "rna_identity"),
                            *(str(fingerprints / f'{r["sample_id"]}.somalier') for r in manifest)],
                           env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
        print("All 14 RNA fingerprints related; evaluate all 91 pair comparisons before interpreting identity.", flush=True)


if __name__ == "__main__":
    main()
