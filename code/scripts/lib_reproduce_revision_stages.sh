#!/usr/bin/env bash
# Called from a verified workspace by reproduce_revision.sh.
set -euo pipefail
Rscript scripts/00_check_revision_environment.R
python3 scripts/verify_replayed_figures.py start
Rscript scripts/verify_retained_figure_inputs.R
Rscript scripts/verify_retained_score_endpoints.R
Rscript scripts/11_scrna_qc_figures.R --revision-replay
AIS_CELL_REPLAY="$(Rscript scripts/rebuild_cell_audits.R)"
python3 scripts/verify_cell_audit_replay.py --root . --run "$AIS_CELL_REPLAY" --report "$AIS_CELL_REPLAY/independent_cell_audit_comparison.json"
Rscript scripts/refit_retained_cell_mixed_models.R
Rscript scripts/refit_specificity_models.R
python3 analysis/scripts/audit_clinical_identity.py
Rscript analysis/bulk_revision_statistics.R
Rscript analysis/bulk_cohort_matched_random_null.R
Rscript analysis/signature_robustness.R
Rscript analysis/signature_pseudobulk_expression.R
python3 scripts/25_audit_deposit_and_clinical_linkage.py
Rscript analysis/verify_deposit_comparisons.R
Rscript analysis/remaining_evidence_sensitivities.R
python3 scripts/verify_remaining_evidence.py
Rscript scripts/refit_full_concordance.R
python3 scripts/verify_fresh_concordance.py
Rscript analysis/bulk_gate_genomewide.R
python3 scripts/gate_candidate_exact.py
python3 scripts/verify_gate_concordance.py --require-generated-normalization
python3 scripts/verify_biology_cells.py
Rscript analysis/scripts/figure_reviewer_corrections.R
Rscript analysis/scripts/figure_real_cell_expression.R --input="$AIS_CELL_REPLAY/analysis/signature_expression_by_library_lineage_before_after.tsv"
Rscript analysis/scripts/figure_rna_pairing.R
Rscript analysis/scripts/figure_rna_failure_diagnostics.R
Rscript analysis/scripts/figure_gate_concordance.R
Rscript scripts/82_additional_validation_figure.R --revision-replay
Rscript scripts/83_gse254300_cs_specificity_figure.R --revision-replay
python3 scripts/verify_replayed_figures.py verify
Rscript analysis/bulk_test_family_inventory.R
python3 scripts/verify_test_family_inventory.py
