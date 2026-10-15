"""Fresh producer outputs must never be replaced by an archived alias."""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from project_paths import ANALYSIS_INPUT_ALIASES, resolve_analysis_input


class AnalysisInputAliasTests(unittest.TestCase):
    legacy = "results/tables/bulk/GSE254300_AIS_paired_edgeR_full.tsv.gz"
    canonical = ANALYSIS_INPUT_ALIASES[legacy]

    def fixture(self, root, name, text):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def test_canonical_read_does_not_restore_duplicate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixed = self.fixture(root, self.canonical, "baseline")
            self.assertEqual(resolve_analysis_input(self.legacy, root), fixed.resolve())
            self.assertFalse((root/self.legacy).exists())

    def test_changed_fresh_output_wins_without_baseline_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixed = self.fixture(root, self.canonical, "baseline")
            fresh = self.fixture(root, self.legacy, "new raw result differs")
            self.assertEqual(resolve_analysis_input(self.legacy, root), fresh.resolve())
            self.assertEqual(fixed.read_text(), "baseline")

    def test_missing_both_sources_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                resolve_analysis_input(self.legacy, Path(tmp))

    def test_path_escape_and_absolute_paths_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ("../escape", "/tmp/escape", "analysis/../escape", "analysis//x", "C:/escape", "analysis\\x"):
                with self.subTest(name=name), self.assertRaises(ValueError):
                    resolve_analysis_input(name, Path(tmp))

    def test_symlink_fresh_file_never_falls_back_to_baseline(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); fixed = self.fixture(root, self.canonical, "baseline")
            fresh = root/self.legacy; fresh.parent.mkdir(parents=True)
            fresh.symlink_to(fixed)
            with self.assertRaisesRegex(ValueError, "Symlink"):
                resolve_analysis_input(self.legacy, root)

    def test_symlink_parent_and_canonical_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root/"outside").mkdir()
            (root/"analysis").symlink_to(root/"outside", target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "Symlink"):
                resolve_analysis_input(self.legacy, root)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript unavailable")
    def test_R_resolver_agrees_on_fresh_baseline_and_security(self):
        # No scientific packages or data are needed for these path checks.
        program = r'''
args <- commandArgs(trailingOnly=TRUE)
source(args[1]); root <- args[2]
legacy <- names(analysis_input_aliases)[1]; canonical <- unname(analysis_input_aliases[legacy])
dir.create(dirname(file.path(root,canonical)), recursive=TRUE)
writeLines("baseline",file.path(root,canonical))
stopifnot(analysis_input(legacy,root)==normalizePath(file.path(root,canonical)), !file.exists(file.path(root,legacy)))
dir.create(dirname(file.path(root,legacy)), recursive=TRUE)
writeLines("different fresh result",file.path(root,legacy))
stopifnot(analysis_input(legacy,root)==normalizePath(file.path(root,legacy)), readLines(file.path(root,canonical))=="baseline")
for (bad in c("../escape", "/tmp/escape", "analysis/../escape", "analysis//x", "C:/escape", "analysis\\x"))
  stopifnot(inherits(try(analysis_input(bad,root),silent=TRUE),"try-error"))
unlink(file.path(root,legacy)); stopifnot(file.symlink(file.path(root,canonical),file.path(root,legacy)))
stopifnot(inherits(try(analysis_input(legacy,root),silent=TRUE),"try-error"))
'''
        with tempfile.TemporaryDirectory() as tmp:
            script = Path(tmp)/"check.R"; script.write_text(program)
            subprocess.run([shutil.which("Rscript"), str(script), str(ROOT/"scripts/analysis_paths.R"), tmp],
                           check=True, capture_output=True, text=True)


if __name__ == "__main__":
    unittest.main()
