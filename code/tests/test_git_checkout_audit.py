"""Git checkout metadata must not weaken the release's scientific file audit."""
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import audit_release as audit


class GitCheckoutAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.parent = Path(self.temp.name)
        self.root = self.parent / 'release'
        self.root.mkdir()
        (self.root / 'README.md').write_text('Fixed release fixture.\n')
        self.write_manifest()

    def tearDown(self):
        self.temp.cleanup()

    def write_manifest(self, extra=()):
        with (self.root / 'SHA256SUMS.tsv').open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['path', 'bytes', 'sha256'])
            for name in ['README.md', *extra]:
                path = self.root / name
                writer.writerow([name, path.stat().st_size,
                                 hashlib.sha256(path.read_bytes()).hexdigest()])

    def git(self, root, *arguments):
        return subprocess.run(['git', '-C', str(root), *arguments], check=True,
                              capture_output=True, text=True)

    def commit_fixture(self):
        self.git(self.root, 'init', '-q')
        self.git(self.root, 'add', 'README.md', 'SHA256SUMS.tsv')
        self.git(self.root, '-c', 'user.name=Release QA', '-c',
                 'user.email=release-qa@example.invalid', 'commit', '-qm', 'Fixture')

    def check(self, root=None):
        # These tests exercise actual Git metadata and file-integrity checks.
        # R source parsing is tested by the complete release audit separately.
        real_run = subprocess.run

        def run(command, *args, **kwargs):
            if command[0] == 'Rscript':
                return subprocess.CompletedProcess(command, 0)
            return real_run(command, *args, **kwargs)

        with patch.object(audit.subprocess, 'run', side_effect=run):
            return audit.audit(self.root if root is None else root)

    def test_actual_git_init_keeps_release_hash_checks(self):
        self.git(self.root, 'init', '-q')
        self.assertEqual(self.check()['manifest_files_verified'], 1)

    def test_actual_local_clone_keeps_release_hash_checks(self):
        self.commit_fixture()
        clone = self.parent / 'clone'
        self.git(self.parent, 'clone', '-q', '--no-hardlinks', str(self.root), str(clone))
        self.assertEqual(self.check(clone)['manifest_files_verified'], 1)

    def test_actual_worktree_gitfile_has_reciprocal_binding(self):
        self.commit_fixture()
        worktree = self.parent / 'worktree'
        self.git(self.root, 'worktree', 'add', '--detach', '-q', str(worktree))
        self.assertTrue((worktree / '.git').is_file())
        self.assertEqual(self.check(worktree)['manifest_files_verified'], 1)

    def test_root_git_symlink_is_never_an_exemption(self):
        other = self.parent / 'other'
        other.mkdir()
        self.git(other, 'init', '-q')
        (self.root / '.git').symlink_to(other / '.git', target_is_directory=True)
        with self.assertRaisesRegex(AssertionError, 'must not be a symlink'):
            self.check()

    def test_invalid_root_git_directory_and_file_are_rejected(self):
        metadata = self.root / '.git'
        metadata.mkdir()
        (metadata / 'hidden.tsv').write_text('unverified\n')
        with self.assertRaisesRegex(AssertionError, 'Invalid root Git metadata'):
            self.check()
        shutil.rmtree(metadata)
        metadata.write_text('gitdir: /not/a/repository\n')
        with self.assertRaisesRegex(AssertionError, 'Invalid root Git metadata'):
            self.check()

    def test_gitfile_cannot_borrow_an_unrelated_git_directory(self):
        other = self.parent / 'other'
        other.mkdir()
        self.git(other, 'init', '-q')
        (self.root / '.git').write_text('gitdir: ' + str(other / '.git') + '\n')
        with self.assertRaisesRegex(AssertionError, 'Invalid root Git metadata'):
            self.check()

    def test_nested_git_directory_including_generated_logs_is_rejected(self):
        for relative in ['nested/.git', 'logs/.git']:
            with self.subTest(relative=relative):
                metadata = self.root / relative
                metadata.mkdir(parents=True)
                with self.assertRaisesRegex(AssertionError, 'Nested Git metadata'):
                    self.check()
                metadata.rmdir()

    def test_nested_gitfile_is_rejected(self):
        folder = self.root / 'logs'
        folder.mkdir()
        (folder / '.git').write_text('gitdir: /not/an/exemption\n')
        with self.assertRaisesRegex(AssertionError, 'Nested Git metadata'):
            self.check()

    def test_git_checkout_does_not_hide_other_unmanifested_files(self):
        self.git(self.root, 'init', '-q')
        (self.root / 'surprise.tsv').write_text('unverified\n')
        with self.assertRaisesRegex(AssertionError, 'Unmanifested release file'):
            self.check()

    def test_symlink_and_broken_symlink_in_generated_paths_are_rejected(self):
        folder = self.root / 'logs'
        folder.mkdir()
        outside = self.parent / 'outside'
        outside.mkdir()
        (outside / 'source.tsv').write_text('unverified\n')
        for target in [outside / 'source.tsv', outside, outside / 'missing']:
            with self.subTest(target=str(target)):
                link = folder / 'link'
                link.symlink_to(target, target_is_directory=target.is_dir())
                with self.assertRaisesRegex(AssertionError, 'Release symlink'):
                    self.check()
                link.unlink()

    def test_manifest_cannot_treat_git_metadata_as_scientific_input(self):
        self.git(self.root, 'init', '-q')
        self.write_manifest(['.git/HEAD'])
        with self.assertRaisesRegex(AssertionError, 'Unsafe/duplicate manifest path'):
            self.check()

    def test_git_checkout_still_rejects_changed_manifest_source(self):
        self.git(self.root, 'init', '-q')
        (self.root / 'README.md').write_text('Different release data.\n')
        with self.assertRaisesRegex(AssertionError, 'Missing/changed release file|SHA-256 mismatch'):
            self.check()


if __name__ == '__main__':
    unittest.main()
