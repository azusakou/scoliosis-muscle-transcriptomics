#!/usr/bin/env python3
"""Check all 16 current-run figure outputs, not PDFs copied as release inputs.

The fixed R producers run between start and verify. Freshness, source exclusion,
file structure, dimensions and hashes are checked without a PDF dependency.
This is an output-production receipt, not a new upstream scientific analysis.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
CONFIG = "config/replayed_figures.json"
MARKER = "logs/figure_replay/started.json"
REPORT = "logs/figure_replay/generated_pdfs.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_path(root, relative):
    name = PurePosixPath(relative)
    if name.is_absolute() or any(p in ("", ".", "..") for p in name.parts) or "\\" in relative:
        raise ValueError("Unsafe replay figure path")
    p = root
    for part in name.parts:
        p = p / part
        if p.is_symlink():
            raise ValueError("Symlink replay figure/producer path: " + relative)
    if root not in p.resolve().parents:
        raise ValueError("Replay figure escapes workspace")
    return p


def configuration(root):
    path = checked_path(root, CONFIG)
    obj = json.loads(path.read_text())
    figures = obj["figures"]
    if len(figures) != 16 or len({f["path"] for f in figures}) != 16 or len({f["id"] for f in figures}) != 16:
        raise ValueError("Expected exactly 16 distinct configured PDF outputs")
    for f in figures:
        if not f["path"].startswith("analysis/figures/") or not f["path"].endswith(".pdf"):
            raise ValueError("Unexpected figure output location")
        checked_path(root, f["path"])
        if not checked_path(root, f["producer"]).is_file():
            raise ValueError("Missing configured figure producer")
    manifest = root / "SHA256SUMS.tsv"
    if manifest.exists():
        with manifest.open(newline="") as handle:
            supplied = {r["path"] for r in csv.DictReader(handle, delimiter="\t")}
        if supplied & {f["path"] for f in figures}:
            raise ValueError("Generated figures must not be copied immutable release inputs")
    return figures, digest(path)


def start(root=ROOT):
    root = root.resolve()
    figures, config_sha = configuration(root)
    obj = {
        "started_ns": time.time_ns(),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "config_sha256": config_sha,
        "producer_sha256": {f["producer"]: digest(checked_path(root, f["producer"])) for f in figures},
        "scope": "Start of current retained plotting run; supplied PDFs are not accepted as generated outputs",
    }
    marker = checked_path(root, MARKER)
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(obj, indent=2) + "\n")
    return obj


def verify(root=ROOT):
    root = root.resolve()
    figures, config_sha = configuration(root)
    marker = checked_path(root, MARKER)
    if not marker.is_file():
        raise ValueError("Missing current figure-production start marker")
    begun = json.loads(marker.read_text())
    if begun["config_sha256"] != config_sha:
        raise ValueError("Figure configuration changed during replay")
    outputs = []
    for f in figures:
        producer = checked_path(root, f["producer"])
        if begun["producer_sha256"].get(f["producer"]) != digest(producer):
            raise ValueError("Figure producer changed during replay")
        path = checked_path(root, f["path"])
        if not path.is_file() or path.stat().st_size < 1024:
            raise ValueError("Missing or empty generated PDF: " + f["path"])
        if path.stat().st_mtime_ns < begun["started_ns"]:
            raise ValueError("Stale/copied PDF predates current producer run: " + f["path"])
        data = path.read_bytes()
        box = re.search(rb"/MediaBox\s*\[\s*([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s*\]", data)
        xref = re.search(rb"startxref\s+(\d+)\s+%%EOF\s*$", data)
        if (not data.startswith(b"%PDF-") or not box or not xref or
                data[int(xref.group(1)):int(xref.group(1))+4] != b"xref" or
                not re.search(rb"/Root\s+\d+\s+0\s+R", data[-4096:]) or
                len(re.findall(rb"/Type\s*/Page\b", data)) != 1):
            raise ValueError("Invalid/incomplete R PDF output: " + f["path"])
        x0, y0, x1, y1 = map(float, box.groups())
        width, height = (x1-x0)*25.4/72, (y1-y0)*25.4/72
        if abs(width-f["width_mm"]) > .5 or abs(height-f["height_mm"]) > .5:
            raise ValueError("Unexpected generated PDF dimensions: " + f["path"])
        outputs.append(dict(f, bytes=len(data), sha256=hashlib.sha256(data).hexdigest(),
                            modified_ns=path.stat().st_mtime_ns, observed_width_mm=width,
                            observed_height_mm=height, producer_sha256=digest(producer)))
    result = {
        "status": "PASS", "generated_PDFs": len(outputs),
        "scope": "All 16 expected PDFs are valid, nonempty current-run outputs, absent from immutable source inputs; no new upstream raw analysis inferred",
        "start_marker": MARKER, "config_sha256": config_sha, "outputs": outputs,
    }
    report = checked_path(root, REPORT)
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["start", "verify"])
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    result = start(args.root) if args.mode == "start" else verify(args.root)
    print(json.dumps({"status": "STARTED", "expected_PDFs": 16, "marker": MARKER}
                     if args.mode == "start" else
                     {"status": result["status"], "generated_PDFs": result["generated_PDFs"], "report": REPORT}, indent=2))
