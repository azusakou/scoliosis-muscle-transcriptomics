#!/usr/bin/env python3
"""Independently check the fresh 32x128 model-refit outputs from this run."""
import hashlib
import json
from pathlib import Path

from verify_edgeR_concordance_null import verify

ROOT = Path(__file__).resolve().parents[1]
run = Path((ROOT / "logs/full_concordance_refit/latest_run.txt").read_text().strip()).resolve()
assert (ROOT / "logs/full_concordance_refit").resolve() in run.parents
receipt = json.loads((run / "model_refit_receipt.json").read_text())
assert receipt["status"] == "MODEL_REFITS_PASSED_PENDING_INDEPENDENT_CHECK"
assert receipt["source_script_SHA256"] == hashlib.sha256((ROOT / "scripts/refit_full_concordance.R").read_bytes()).hexdigest()
assert receipt["source_bindings_SHA256"] == hashlib.sha256((ROOT / "analysis/bulk_gate_inputs/source_bindings.tsv").read_bytes()).hexdigest()
result = verify(run)
result["scope"] = "Fresh count-model refits and independent 4096-grid arithmetic; no raw-read realignment or anatomical authentication"
(run / "independent_concordance.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
