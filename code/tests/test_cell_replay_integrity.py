"""Source-identity and provenance guards for independent fresh cell-audit checks."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import verify_cell_audit_replay as replay


class CellReplayIntegrityTests(unittest.TestCase):
    def fixture(self, directory):
        root = Path(directory) / "release"; root.mkdir()
        run = Path(directory) / "fresh"; run.mkdir()
        script = root / "scripts/rebuild_cell_audits.R"
        script.parent.mkdir(); script.write_text("# fixed rebuilder\n")
        for parent in (root, run):
            metadata = parent / "results/tables/annotation/annotated_cell_metadata.tsv.gz"
            metadata.parent.mkdir(parents=True); metadata.write_bytes(b"fixed retained metadata")
        record = dict(rebuild_script_sha256=replay.reproduction.digest(script),
                      metadata_sha256=replay.reproduction.digest(root / "results/tables/annotation/annotated_cell_metadata.tsv.gz"),
                      input_files=[])
        return root, run, record

    def save(self, run, record):
        (run / "cell_audit_replay.json").write_text(json.dumps(record))

    def test_changed_rebuilder_cannot_authenticate_an_earlier_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root, run, record = self.fixture(directory); self.save(run, record)
            (root / "scripts/rebuild_cell_audits.R").write_text("# changed producer\n")
            with self.assertRaisesRegex(ValueError, "different rebuilder version"):
                replay.verify(root, run)

    def test_changed_fresh_cell_metadata_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root, run, record = self.fixture(directory); self.save(run, record)
            (run / "results/tables/annotation/annotated_cell_metadata.tsv.gz").write_bytes(b"wrong labels or totals")
            with self.assertRaisesRegex(ValueError, "metadata must match"):
                replay.verify(root, run)

    def test_changed_sparse_input_cannot_reuse_original_input_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root, run, record = self.fixture(directory)
            source = root / "data/retained/cell_target_counts.tsv.gz"
            source.parent.mkdir(parents=True); source.write_bytes(b"original counts")
            record["input_files"] = [dict(path=source.name, bytes=source.stat().st_size,
                                          sha256=replay.reproduction.digest(source))]
            self.save(run, record); source.write_bytes(b"modified counts")
            with self.assertRaisesRegex(ValueError, "input differs"):
                replay.verify(root, run)

    def manifests(self, directory, fresh_source="data/retained/cell_target_counts.tsv.gz", feature_count=38224):
        roots = (Path(directory) / "reference", Path(directory) / "fresh")
        fields = ("assay", "n_genes", "n_cells", "source_object")
        for root, source, count in ((roots[0], "data/processed/scrna/scrna_annotated.rds", 38224),
                                   (roots[1], fresh_source, feature_count)):
            path = root / "analysis/signature_expression_assay_manifest.tsv"
            path.parent.mkdir(parents=True)
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(fields)
                writer.writerows((assay, count, 14846, source) for assay in ("RNA", "SoupX"))
        return roots

    def test_sparse_replay_cannot_claim_it_loaded_the_historical_rds(self):
        with tempfile.TemporaryDirectory() as directory:
            root, run = self.manifests(directory, "data/processed/scrna/scrna_annotated.rds")
            with self.assertRaisesRegex(ValueError, "declare its sparse-count source"):
                replay.assay_manifest_comparison(root, run)

    def test_full_assay_feature_count_remains_exact(self):
        with tempfile.TemporaryDirectory() as directory:
            root, run = self.manifests(directory, feature_count=38225)
            with self.assertRaisesRegex(ValueError, "Changed integer/count value"):
                replay.assay_manifest_comparison(root, run)


if __name__ == "__main__":
    unittest.main()
