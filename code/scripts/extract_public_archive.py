#!/usr/bin/env python3
"""Extract regular public-data files only, preserving safe relative member paths."""
from __future__ import annotations

import argparse
import os
import shutil
import stat
import tarfile
import tempfile
import zipfile
from pathlib import Path, PurePosixPath


def destination_path(root: Path, name: str) -> Path:
    member = PurePosixPath(name)
    if not name or "\\" in name or ":" in name or member.is_absolute() or ".." in member.parts:
        raise ValueError(f"Unsafe archive member: {name!r}")
    destination = (root / Path(*member.parts)).resolve()
    if destination == root.resolve() or root.resolve() not in destination.parents:
        raise ValueError(f"Archive member escapes destination: {name!r}")
    return destination


def _copy_file(source, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as output:
            temporary = Path(output.name)
            shutil.copyfileobj(source, output, length=1024 * 1024)
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def extract_archive(archive: Path, destination: Path) -> int:
    destination.mkdir(parents=True, exist_ok=True)
    count = 0
    seen = set()
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as handle:
            members = []
            for item in handle.infolist():
                if PurePosixPath(item.filename).name.startswith("._"):
                    continue
                target = destination_path(destination, item.filename.rstrip("/"))
                if stat.S_IFMT(item.external_attr >> 16) not in (0, stat.S_IFREG, stat.S_IFDIR):
                    raise ValueError(f"Archive links are not accepted: {item.filename}")
                if target in seen:
                    raise ValueError(f"Duplicate archive destination: {item.filename}")
                seen.add(target)
                members.append((item, target))
            for item, target in members:
                if item.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                else:
                    with handle.open(item) as source:
                        _copy_file(source, target)
                    count += 1
    elif tarfile.is_tarfile(archive):
        with tarfile.open(archive) as handle:
            members = []
            for item in handle.getmembers():
                if PurePosixPath(item.name).name.startswith("._"):
                    continue
                # Tar root directory entries carry no data.
                if item.isdir() and item.name.rstrip("/") in (".", ""):
                    continue
                target = destination_path(destination, item.name.rstrip("/"))
                if not (item.isdir() or item.isfile()):
                    raise ValueError(f"Archive links/devices are not accepted: {item.name}")
                if target in seen:
                    raise ValueError(f"Duplicate archive destination: {item.name}")
                seen.add(target)
                members.append((item, target))
            for item, target in members:
                if item.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                else:
                    source = handle.extractfile(item)
                    if source is None:
                        raise ValueError(f"Cannot read archive member: {item.name}")
                    with source:
                        _copy_file(source, target)
                    count += 1
    else:
        raise ValueError(f"Not a supported ZIP/TAR data archive: {archive}")
    return count


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    print(f"Extracted {extract_archive(args.archive, args.destination)} regular files")
