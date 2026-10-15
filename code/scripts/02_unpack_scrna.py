#!/usr/bin/env python3
"""Unpack public processed/raw feature-barcode bundles by detected content."""

import argparse
import pathlib
from extract_public_archive import extract_archive
from project_paths import raw_data_root


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=pathlib.Path, default=raw_data_root())
    args = parser.parse_args()
    source = args.raw_root / "scrna"
    destination = source / "matrices"
    destination.mkdir(parents=True, exist_ok=True)
    for archive in sorted(source.glob("*cellranger*")):
        if not archive.is_file():
            continue
        count = extract_archive(archive, destination)
        print(f"OK\t{archive}\t{count} files")


if __name__ == "__main__":
    main()

