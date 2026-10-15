#!/usr/bin/env python3
from __future__ import annotations

from project_paths import raw_data_root, executable, curl_retry_options

import csv
import hashlib
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "metadata/audit/gse311180_requant_manifest.tsv"
RAW_ROOT = raw_data_root()
OUT = RAW_ROOT / "bulk_gse311180" / "sra"
CURL = executable("curl", "AIS_CURL")
CURL_RETRIES = curl_retry_options(CURL)
LOG = ROOT / "logs/gse311180_sra_download.tsv"
if not RAW_ROOT.exists():
    raise SystemExit(f"Configured AIS_RAW_ROOT does not exist: {RAW_ROOT}")
OUT.mkdir(parents=True, exist_ok=True)
LOG.parent.mkdir(parents=True, exist_ok=True)

with MANIFEST.open() as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def locate(run: str) -> Path | None:
    candidates = [OUT / run / f"{run}.sra", OUT / run / run, OUT / f"{run}.sra"]
    return next((p for p in candidates if p.exists()), None)


def download(row: dict[str, str]) -> dict[str, str | int]:
    run = row["sra_run"]
    if not RAW_ROOT.exists():
        return {"sra_run": run, "status": "EXTERNAL_ROOT_UNAVAILABLE", "path": "",
                "bytes": 0, "md5": "", "expected_md5": row["sra_md5"],
                "message": f"AIS_RAW_ROOT unavailable during transfer: {RAW_ROOT}"}
    expected_size = int(row["sra_size_bytes"])
    run_dir = OUT / run
    recovery_message = lambda: "quarantined=" + ",".join(
        sorted(p.name for p in run_dir.glob(f"{run}.sra.invalid_*")))
    existing = locate(run)
    if existing is not None:
        observed_size = existing.stat().st_size
        observed_md5 = md5(existing) if observed_size == expected_size else ""
        if observed_size == expected_size and observed_md5 == row["sra_md5"]:
            return {"sra_run": run, "status": "PASS", "path": str(existing),
                    "bytes": observed_size, "md5": observed_md5,
                    "expected_md5": row["sra_md5"], "message": recovery_message()}
        reason = f"size_{observed_size}" if observed_size != expected_size else f"md5_{observed_md5[:12]}"
        existing.replace(existing.with_name(f"{existing.name}.invalid_{reason}"))
        existing = None
    if existing is None:
        run_dir.mkdir(parents=True, exist_ok=True)
        partial = run_dir / f"{run}.sra.tmp"
        if partial.exists() and partial.stat().st_size > expected_size:
            # A stale duplicate downloader can only make a resumable object
            # larger than the repository object. Quarantine it rather than
            # deleting it, then restart this one accession from byte zero.
            observed_size = partial.stat().st_size
            quarantine = run_dir / f"{run}.sra.invalid_oversized_{observed_size}"
            partial.replace(quarantine)
        proc = subprocess.run(
            [CURL, "--fail", "--location", "--silent", "--show-error",
             *CURL_RETRIES, "--continue-at", "-",
             "--output", str(partial), row["sra_url"]],
            text=True, capture_output=True,
        )
        if proc.returncode != 0:
            return {"sra_run": run, "status": "DOWNLOAD_FAILED", "path": "",
                    "bytes": 0, "md5": "", "expected_md5": row["sra_md5"],
                    "message": (proc.stderr or proc.stdout)[-1000:]}
        if partial.stat().st_size == expected_size:
            partial.replace(run_dir / f"{run}.sra")
        else:
            return {"sra_run": run, "status": "SIZE_MISMATCH_AFTER_DOWNLOAD",
                    "path": str(partial), "bytes": partial.stat().st_size,
                    "md5": "", "expected_md5": row["sra_md5"],
                    "message": f"expected {expected_size} bytes"}
    path = locate(run)
    if path is None:
        return {"sra_run": run, "status": "FILE_NOT_FOUND_AFTER_PREFETCH", "path": "",
                "bytes": 0, "md5": "", "expected_md5": row["sra_md5"], "message": ""}
    observed = md5(path)
    ok = path.stat().st_size == expected_size and observed == row["sra_md5"]
    return {"sra_run": run, "status": "PASS" if ok else "CHECKSUM_OR_SIZE_FAIL",
            "path": str(path), "bytes": path.stat().st_size, "md5": observed,
            "expected_md5": row["sra_md5"], "message": recovery_message()}


results: list[dict[str, str | int]] = []
# Four parallel transfers saturate the available connection without making the
# external volume contend with simultaneous fasterq/Salmon I/O. Interrupted
# `.sra.tmp` files are resumed by curl and accepted only after size and MD5
# validation.
with ThreadPoolExecutor(max_workers=4) as pool:
    futures = {pool.submit(download, row): row["sra_run"] for row in rows}
    for future in as_completed(futures):
        result = future.result()
        results.append(result)
        print(f"{result['sra_run']}\t{result['status']}", flush=True)

results.sort(key=lambda x: str(x["sra_run"]))
with LOG.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

failed = [r for r in results if r["status"] != "PASS"]
if failed:
    raise SystemExit(f"{len(failed)} SRA downloads failed validation; see {LOG}")
print(f"All {len(results)} SRA objects passed size and MD5 validation")
