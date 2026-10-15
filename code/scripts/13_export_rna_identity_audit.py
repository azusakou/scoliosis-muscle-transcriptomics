#!/usr/bin/env python3
"""Export completed full-read RNA pairing/QC aggregates, never genotype sketches.

Run full source validation (08) and all-sample fingerprint validation (12) first.
This checks their provenance, uniform recipes, completed quantifications, native
fingerprint hashes, coverage and complete read counts before publishing tables.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import tempfile
from input_integrity import digest, validate_quantification
from project_paths import raw_data_root
from verify_rna_identity import assess, table, verify_mapping_qc

ROOT = Path(__file__).resolve().parents[1]


def write_table(destination, name, rows):
    if not rows or any(set(row) != set(rows[0]) for row in rows):
        raise ValueError("Empty or inconsistent aggregate schema: " + name)
    with tempfile.NamedTemporaryFile(mode="w", newline="", encoding="utf-8",
                                     dir=destination, delete=False) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, destination / name)


def export(identity, raw, quant, destination, manifest_path, criteria_path):
    manifest = table(manifest_path)
    criteria = json.loads(criteria_path.read_text())
    pair_source = identity / "rna_identity.pairs.tsv"
    sample_source = identity / "rna_identity.samples.tsv"
    binding_path = identity / "reference_artifact_binding.json"
    reference_binding = json.loads(binding_path.read_text())
    pairs, individuals = table(pair_source), table(sample_source)
    assessment = assess(pairs, manifest, criteria)
    if (len(individuals) != 14 or len({r["sample_id"] for r in individuals}) != 14 or
            {r["sample_id"] for r in individuals} != {r["sample_id"] for r in manifest}):
        raise ValueError("Incomplete/duplicate individual fingerprint QC")
    individuals = {r["sample_id"]: r for r in individuals}
    canonical = []
    for row in pairs:
        value = dict(row)
        value["somalier_reported_relatedness"] = value["relatedness"]
        value["relatedness"] = 2 * (int(value["shared_hets"]) - 2 * int(value["ibs0"])) / max(1, int(value["hets_ab"]))
        canonical.append(value)
    qc, provenance, recipes = [], [], set()
    for row in manifest:
        sample = row["sample_id"]
        fqpath = raw / "bulk_gse311180/fastq" / f"{sample}.provenance.json"
        fq = json.loads(fqpath.read_text())
        qpath = quant / sample
        q = validate_quantification(qpath, expected_read_pairs=fq["outputs"]["read_pairs"])
        qstamp = json.loads((quant / f"{sample}.provenance.json").read_text())
        if qstamp["inputs"]["fastq"] != fq["outputs"] or qstamp["outputs"] != q:
            raise ValueError("Quantification provenance mismatch: " + sample)
        meta = json.loads((qpath / "aux_info/meta_info.json").read_text())
        stamp = json.loads((identity / f"{sample}.provenance.json").read_text())
        if stamp["inputs"]["fastq_provenance_sha256"] != digest(fqpath):
            raise ValueError("Fingerprint source mismatch: " + sample)
        recipes.add(json.dumps(stamp["inputs"]["recipe"], sort_keys=True))
        if stamp["inputs"]["recipe"]["script_sha256"] != reference_binding["binding"]["scientific_producer_sha256"]:
            raise ValueError("Fingerprint/reference scientific-producer mismatch: " + sample)
        fingerprint = identity / "fingerprints" / f"{sample}.somalier"
        if digest(fingerprint) != stamp["outputs"]["fingerprint_sha256"]:
            raise ValueError("Completed fingerprint hash mismatch: " + sample)
        logpath = identity / sample / "Log.final.out"
        if digest(logpath) != stamp["outputs"]["STAR_log_sha256"]:
            raise ValueError("Completed STAR log mismatch: " + sample)
        star = {line.split("|", 1)[0].strip(): line.split("|", 1)[1].strip()
                for line in logpath.read_text().splitlines() if "|" in line}
        n = int(fq["outputs"]["read_pairs"])
        if int(star["Number of input reads"]) != n or int(meta["num_processed"]) != n:
            raise ValueError("Complete input read count mismatch: " + sample)
        called = sum(int(individuals[sample][k]) for k in ("n_hom_ref", "n_het", "n_hom_alt"))
        qc.append(dict(sample_id=sample, patient_id=row["patient_id"], side=row["side"], sra_run=row["sra_run"],
            validated_read_pairs=n, salmon_processed_pairs=int(meta["num_processed"]),
            salmon_mapped_pairs=int(meta["num_mapped"]), salmon_percent_mapped=meta["percent_mapped"],
            STAR_input_pairs=int(star["Number of input reads"]),
            STAR_uniquely_mapped_pairs=int(star["Uniquely mapped reads number"]),
            STAR_uniquely_mapped_percent=star["Uniquely mapped reads %"].rstrip("%"),
            STAR_multimapped_pairs=int(star["Number of reads mapped to multiple loci"]),
            somalier_called_sites=called, somalier_unknown_sites=int(individuals[sample]["n_unknown"]),
            somalier_genotype_depth_mean=individuals[sample]["gt_depth_mean"],
            somalier_n_het=int(individuals[sample]["n_het"])))
        provenance.append(dict(sample_id=sample, sra_run=row["sra_run"],
            official_sra_bytes=row["sra_size_bytes"], official_sra_md5=row["sra_md5"],
            mate1_sha256=fq["outputs"]["mate1_sha256"], mate2_sha256=fq["outputs"]["mate2_sha256"],
            FASTQ_provenance_sha256=digest(fqpath), quant_sha256=q["quant_sha256"],
            quant_metadata_sha256=q["meta_sha256"], RNA_BAM_sha256=stamp["outputs"]["bam_sha256"],
            RNA_fingerprint_sha256=stamp["outputs"]["fingerprint_sha256"],
            STAR_log_sha256=stamp["outputs"]["STAR_log_sha256"],
            RNA_recipe_script_sha256=stamp["inputs"]["recipe"]["script_sha256"]))
    if len(recipes) != 1:
        raise ValueError("Mixed RNA alignment/fingerprint recipes; rerun every library uniformly")
    assessment["full_read_mapping_QC_checks"] = verify_mapping_qc(qc, manifest, canonical)
    assessment["uniform_RNA_recipe"] = json.loads(next(iter(recipes)))
    assessment["source_sha256"] = {"native_all_pairs": digest(pair_source), "native_sample_QC": digest(sample_source),
                                   "criteria": digest(criteria_path), "manifest": digest(manifest_path),
                                   "export_script": digest(Path(__file__)),
                                   "reference_artifact_binding": digest(binding_path)}
    destination.mkdir(parents=True, exist_ok=True)
    for name, rows in [("audit_gse311180_rna_all_pairs.tsv", canonical),
                       ("audit_gse311180_rna_expected_pairs.tsv", assessment["expected_pairs"]),
                       ("audit_gse311180_raw_RNA_mapping_QC.tsv", qc),
                       ("audit_gse311180_completed_raw_source_sha256.tsv", provenance)]:
        write_table(destination, name, rows)
    return assessment


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--identity-output", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, default=raw_data_root())
    parser.add_argument("--quant-root", type=Path)
    parser.add_argument("--destination", type=Path, default=ROOT / "analysis")
    parser.add_argument("--manifest", type=Path, default=ROOT / "metadata/audit/gse311180_requant_manifest.tsv")
    parser.add_argument("--criteria", type=Path, default=ROOT / "metadata/audit/rna_identity_consistency_criteria.json")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    quant = args.quant_root or args.raw_root.parent / "processed/gse311180/salmon"
    result = export(args.identity_output, args.raw_root, quant, args.destination, args.manifest, args.criteria)
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(encoded)
    print(encoded)
    print("COMPLETE_RAW_AND_RNA_AGGREGATE_EXPORT_PASS")


if __name__ == "__main__":
    main()
