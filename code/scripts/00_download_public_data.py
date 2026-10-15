#!/usr/bin/env python3
"""Download official metadata and public matrices with auditable URLs.

The downloader is conservative: existing validated files are not overwritten
unless --force is supplied. All data are written below --raw-root.
"""

from __future__ import annotations

from project_paths import raw_data_root
from input_integrity import validate_download

import argparse
import json
import os
import pathlib
import shutil
import sys
import tempfile
import urllib.request
import zipfile


DEFAULT_ROOT = raw_data_root()


def download(url: str, destination: pathlib.Path, force: bool = False) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not force:
        validate_download(destination)
        print(f"SKIP_VALIDATED\t{destination}")
        return
    request = urllib.request.Request(url, headers={"User-Agent": "ais-paraspinal-ecology/1.0"})
    temporary = None
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as handle:
                temporary = pathlib.Path(handle.name)
                shutil.copyfileobj(response, handle)
        validate_download(temporary, destination.name)
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    print(f"OK\t{destination}\t{destination.stat().st_size}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=pathlib.Path, default=DEFAULT_ROOT)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()

    root = args.raw_root
    metadata = {
        "metadata/PRJNA722100_bioproject_esearch.json": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=bioproject&term=PRJNA722100%5BProject%20Accession%5D&retmode=json",
        "metadata/PRJNA722100_bioproject.xml": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=bioproject&id=722100&retmode=xml",
        "metadata/PRJNA722100_sra_esearch.json": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA722100%5BBioProject%5D&retmax=1000&retmode=json",
        "metadata/PRJNA722100_sra_experiments.xml": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=14080942,14080941,14080940,14080939&retmode=xml",
        "metadata/GSE254300_family.soft": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254300&targ=all&view=full&form=text",
        "metadata/GSE254300_supplementary_index.html": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254300/suppl/",
        "metadata/relevant_pubmed_records.xml": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=38440359,41836551,37185898,42244182,30241458,35137485,37401090&retmode=xml",
        "metadata/PMC10130095_original_paper_bioc.json": "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC10130095/unicode",
        "metadata/PMC10910612_bulk_paper_bioc.json": "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC10910612/unicode",
        "bulk_gse311180/GSE311180_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311180/soft/GSE311180_family.soft.gz",
        "bulk_gse311180/PRJNA1368770_sra_metadata.xml": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=41673893,41673892,41673891,41673890,41673889,41673888,41673887,41673886,41673885,41673884,41673883,41673882,41673881,41673880,41673879,41673878,41673877,41673876,41673875,41673874,41673873&retmode=xml",
        "bulk_gse311180/mmc2.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S2666247726000710-mmc2.xlsx",
        "bulk_gse311180/PMC13316586.xml": "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC13316586/unicode",
        "bulk_gse324067/GSE324067_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE324nnn/GSE324067/soft/GSE324067_family.soft.gz",
        "bulk_gse324067/GSE324067_family.xml.tgz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE324nnn/GSE324067/miniml/GSE324067_family.xml.tgz",
        "bulk_gse324067/PRJNA1433384_sra_metadata.xml": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=43284803,43284802,43284801,43284800,43284799,43284798,43284797,43284796,43284795,43284794,43284793,43284792,43284791,43284790,43284789,43284788,43284787,43284786&retmode=xml",
        "reference_gse17503/GSE17503_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17503/soft/GSE17503_family.soft.gz",
        "literature_tent5a_2022/PMC8891553_fullTextXML.xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8891553/fullTextXML",
    }
    data_files = {
            "bulk/GSE254300_readcount.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254300/suppl/GSE254300_readcount.xlsx",
            "bulk/GSE254300_fpkm.xlsx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254300/suppl/GSE254300_fpkm.xlsx",
            "scrna/H01045_concave_cellranger.tar.gz": "https://drive.usercontent.google.com/download?id=1IkX-gLNTf5Tg_0GQ13mKDJh99xk5p9rX&export=download&confirm=t",
            "scrna/H01045_convex_cellranger.tar.gz": "https://drive.usercontent.google.com/download?id=1M4ueXhpUxCef1t-rpwNhTBhYe-73-OnW&export=download&confirm=t",
            "scrna/H66206_concave_cellranger.tar.gz": "https://drive.usercontent.google.com/download?id=1D3utrtoZ3SaEBe4UQwdupJyxG_Qrq3dI&export=download&confirm=t",
            "scrna/H66206_convex_cellranger.tar.gz": "https://drive.usercontent.google.com/download?id=16t7GHJ-CfSaaZ-6SLvf7K9zIqFnbCR1v&export=download&confirm=t",
            "bulk_gse311180/GSE311180_TPM_allSamples.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311180/suppl/GSE311180_TPM_allSamples.txt.gz",
            "bulk_gse311180/mmc3.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S2666247726000710-mmc3.xlsx",
            "bulk_gse311180/mmc5.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S2666247726000710-mmc5.xlsx",
            "bulk_gse311180/mmc6.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S2666247726000710-mmc6.xlsx",
            "bulk_gse311180/mmc15.xlsx": "https://ars.els-cdn.com/content/image/1-s2.0-S2666247726000710-mmc15.xlsx",
            "bulk_gse311180/mmc18.zip": "https://ars.els-cdn.com/content/image/1-s2.0-S2666247726000710-mmc18.zip",
            "bulk_gse324067/GSE324067_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE324nnn/GSE324067/suppl/GSE324067_RAW.tar",
            "reference_gse17503/GSE17503_non-normalized_data.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17503/suppl/GSE17503_non-normalized_data.txt.gz",
            "reference_gse17503/GSE17503_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17503/matrix/GSE17503_series_matrix.txt.gz",
            "literature_h19_2018/10020_2018_49_MOESM1_ESM.docx": "https://static-content.springer.com/esm/art%3A10.1186%2Fs10020-018-0049-y/MediaObjects/10020_2018_49_MOESM1_ESM.docx",
            "literature_tent5a_2022/PMC8891553_supplementaryFiles.zip": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8891553/supplementaryFiles",
        }
    declared_files = dict(metadata)
    declared_files.update(data_files)
    files = metadata if args.metadata_only else declared_files

    for relative, url in files.items():
        destination = root / relative
        download(url, destination, args.force)

    # Europe PMC returns a bundle containing article figures and supplements.
    # Extract only the expected TENT5A DOCX with a fixed member name.
    tent5a_bundle = root / "literature_tent5a_2022/PMC8891553_supplementaryFiles.zip"
    if tent5a_bundle.exists():
        with zipfile.ZipFile(tent5a_bundle) as archive:
            member = "CPR-55-e13183-s002.docx"
            if member not in archive.namelist():
                raise RuntimeError(f"Expected {member} is absent from {tent5a_bundle}")
            destination = root / "literature_tent5a_2022" / member
            if args.force or not destination.exists() or destination.stat().st_size == 0:
                destination.write_bytes(archive.read(member))
                print(f"OK\t{destination}\t{destination.stat().st_size}")
    provenance = []
    for relative, url in declared_files.items():
        destination = root / relative
        provenance.append({
            "path": relative,
            "url": url,
            "bytes": destination.stat().st_size if destination.is_file() else None,
            "local_status": "present" if destination.is_file() else "not_downloaded",
        })
    project_manifest = pathlib.Path(__file__).resolve().parents[1] / "metadata/audit/download_manifest.json"
    project_manifest.parent.mkdir(parents=True, exist_ok=True)
    project_manifest.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    external_manifest = root / "metadata/download_manifest.json"
    try:
        external_manifest.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    except PermissionError:
        print(f"WARN\tCannot mirror manifest to read-only external root: {external_manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
