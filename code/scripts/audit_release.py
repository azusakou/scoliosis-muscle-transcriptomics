#!/usr/bin/env python3
"""Check a prepared public release: completeness, hashes, source syntax and exclusions."""
from __future__ import annotations
import argparse
import ast
import csv
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
from input_integrity import digest


def root_git_metadata(root: Path) -> str | None:
    """Recognize only a real checkout's root metadata, never scientific inputs."""
    metadata = root / '.git'
    if metadata.is_symlink():
        raise AssertionError('Root Git metadata must not be a symlink')
    if not metadata.exists():
        return None
    if not metadata.is_dir() and not metadata.is_file():
        raise AssertionError('Invalid root Git metadata')
    # Do not let inherited Git overrides make another repository authenticate
    # this directory. rev-parse does not execute hooks or access a remote.
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith('GIT_')}
    try:
        result = subprocess.run(
            ['git', '-C', str(root), 'rev-parse', '--absolute-git-dir', '--show-toplevel'],
            check=True, capture_output=True, text=True, env=environment)
        gitdir_text, worktree_text = result.stdout.strip().splitlines()
        gitdir = Path(gitdir_text).resolve()
        if Path(worktree_text).resolve() != root:
            raise ValueError('Git worktree does not match release root')
        if metadata.is_dir():
            if gitdir != metadata.resolve():
                raise ValueError('Git directory does not match root metadata')
            return 'directory'
        # Worktree .git files intentionally point to the parent repository's
        # administrative directory. Require its reciprocal worktree binding;
        # an arbitrary pointer to an unrelated repository is not an exemption.
        pointer = metadata.read_text().strip()
        if not pointer.startswith('gitdir: ') or '\n' in pointer:
            raise ValueError('Malformed worktree Git pointer')
        declared = Path(pointer[len('gitdir: '):])
        if not declared.is_absolute():
            declared = root / declared
        backlink = gitdir / 'gitdir'
        if declared.resolve() != gitdir or not backlink.is_file() or backlink.is_symlink():
            raise ValueError('Missing reciprocal worktree binding')
        if Path(backlink.read_text().strip()).resolve() != metadata:
            raise ValueError('Worktree backlink does not match release root')
        return 'file'
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        raise AssertionError('Invalid root Git metadata') from error


def release_files(root: Path, git_metadata: str | None):
    """Walk without following links; root checkout metadata is the sole exemption."""
    for directory, subdirectories, filenames in os.walk(root, followlinks=False):
        parent = Path(directory)
        for name in list(subdirectories):
            path = parent / name
            if path.is_symlink():
                raise AssertionError(f'Release symlink is not allowed: {path.relative_to(root)}')
            if name == '.git':
                if parent == root and git_metadata == 'directory':
                    subdirectories.remove(name)
                else:
                    raise AssertionError(f'Nested Git metadata is not allowed: {path.relative_to(root)}')
        for name in filenames:
            path = parent / name
            if path.is_symlink():
                raise AssertionError(f'Release symlink is not allowed: {path.relative_to(root)}')
            if name == '.git':
                if parent == root and git_metadata == 'file':
                    continue
                raise AssertionError(f'Nested Git metadata is not allowed: {path.relative_to(root)}')
            if path.is_file():
                yield path


def audit(root: Path, check_manifest: bool = True) -> dict:
    root = root.resolve()
    git_metadata = root_git_metadata(root)
    files = list(release_files(root, git_metadata))
    verified = 0
    if check_manifest:
        with (root / 'SHA256SUMS.tsv').open() as handle:
            records = list(csv.DictReader(handle, delimiter='\t'))
        seen = set()
        for record in records:
            name = record['path']; relative = PurePosixPath(name)
            if relative.is_absolute() or '..' in relative.parts or '.git' in relative.parts or '\\' in name or name in seen:
                raise AssertionError(f'Unsafe/duplicate manifest path: {name}')
            seen.add(name)
            path = root / relative
            if root not in path.resolve().parents or not path.is_file() or path.is_symlink() or path.stat().st_size != int(record['bytes']):
                raise AssertionError(f'Missing/changed release file: {name}')
            if digest(path) != record['sha256']:
                raise AssertionError(f'SHA-256 mismatch: {name}')
            verified += 1
        # Initial contents must be complete; generated validation/reproduction outputs are permitted.
        for path in files:
            name = str(path.relative_to(root))
            generated = name.startswith(('logs/', 'data/processed/', 'analysis/figures/')) or '__pycache__' in path.parts
            if name not in seen and name != 'SHA256SUMS.tsv' and not generated:
                raise AssertionError(f'Unmanifested release file: {name}')
    checked = {'python': 0, 'R': 0, 'shell': 0}
    rfiles = []
    for path in sorted(files):
        if any(part in {'__pycache__', 'logs'} for part in path.relative_to(root).parts):
            continue
        relative = str(path.relative_to(root))
        if path.suffix in {'.pyc', '.pem', '.key'} or 'reply_reviewer' in relative or relative.startswith('bmc/comments/'):
            raise AssertionError(f'Excluded/private artifact found in public package: {relative}')
        if path.suffix in {'.py', '.R', '.sh'}:
            code = path.read_text()
            if re.search(r'/(Volumes|Users)/', code):
                raise AssertionError(f'Personal executable path found: {relative}')
            if re.search(r'(?:gh[pousr]_[A-Za-z0-9]{30,}|sk-[A-Za-z0-9_-]{35,}|-----BEGIN [A-Z ]*PRIVATE KEY-----)', code):
                raise AssertionError(f'Potential credential found: {relative}')
            if path.suffix == '.py':
                ast.parse(code, filename=relative); checked['python'] += 1
            elif path.suffix == '.sh':
                subprocess.run(['bash', '-n', str(path)], check=True); checked['shell'] += 1
            else:
                rfiles.append(str(path)); checked['R'] += 1
    subprocess.run(['Rscript', '-e', 'for (p in commandArgs(TRUE)) parse(file=p)', *rfiles], check=True,
                   stdout=subprocess.DEVNULL)
    return {'status': 'PASS', 'manifest_files_verified': verified, 'syntax_files': checked,
            'scope': 'Prepared release file/syntax checks; not raw-data scientific validation'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--no-manifest', action='store_true')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    result = audit(args.root, not args.no_manifest)
    text = json.dumps(result, indent=2) + '\n'
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True); args.report.write_text(text)
    print(text)
