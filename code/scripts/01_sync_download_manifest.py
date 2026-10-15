#!/usr/bin/env python3
from __future__ import annotations

from project_paths import raw_data_root

import json
import os
from pathlib import Path

project = Path(__file__).resolve().parents[1]
raw_root = raw_data_root()
source = raw_root / "metadata" / "download_manifest.json"
target = project / "metadata" / "audit" / "download_manifest.json"
if not source.is_file():
    if target.is_file():
        print(f"External manifest absent; retained project manifest {target}")
        raise SystemExit(0)
    raise SystemExit(f"Both external and project download manifests are missing: {source}")
records = json.loads(source.read_text(encoding="utf-8"))
if not isinstance(records, list) or not records:
    raise SystemExit(f"Invalid or empty download manifest: {source}")
for record in records:
    if not all(key in record for key in ("path", "url", "bytes")):
        raise SystemExit("Download manifest record lacks path/url/bytes")
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
print(f"Synced {len(records)} provenance records to {target}")
