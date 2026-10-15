#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p logs/release_validation
python3 scripts/audit_release.py --no-manifest --report logs/release_validation/file_syntax_audit.json
python3 -m unittest discover -s tests -v > logs/release_validation/python_tests.log 2>&1
Rscript tests/test_common.R > logs/release_validation/r_tests.log 2>&1
Rscript scripts/verify_retained_models.R > logs/release_validation/retained_model_refits.log 2>&1
python3 scripts/verify_results.py --report logs/release_validation/independent_results.json
python3 scripts/verify_cell_expression.py --report logs/release_validation/independent_cell_expression.json
python3 scripts/verify_edgeR_concordance_null.py --report logs/release_validation/independent_edgeR_null.json
python3 scripts/verify_gate_concordance.py --report logs/release_validation/independent_gate_concordance.json
python3 scripts/verify_biology_cells.py --report logs/release_validation/independent_biological_cells.json
python3 scripts/verify_rna_identity.py --pairs analysis/audit_gse311180_rna_all_pairs.tsv --manifest metadata/audit/gse311180_requant_manifest.tsv --criteria metadata/audit/rna_identity_consistency_criteria.json --expected-pairs analysis/bulk_gate_inputs/RNA_expected_pairs.tsv --mapping-qc analysis/audit_gse311180_raw_RNA_mapping_QC.tsv --output logs/release_validation/independent_RNA_pairing.json
python3 scripts/verify_rna_failure_diagnostics.py --report logs/release_validation/independent_RNA_failure_diagnostics.json
python3 scripts/verify_remaining_evidence.py > logs/release_validation/independent_remaining_evidence.json
python3 scripts/verify_test_family_inventory.py > logs/release_validation/independent_test_inventory.json
Rscript analysis/verify_deposit_comparisons.R > logs/release_validation/independent_deposit_comparisons.log 2>&1
printf 'Release checks passed. Inspect logs/release_validation for the independently checked scope.\n'
