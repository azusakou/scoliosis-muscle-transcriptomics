#!/usr/bin/env Rscript
source("scripts/analysis_paths.R")
# Inventory of retained original analysis families. This counts existing tests;
# it does not pool their P values or impose a retrospective omnibus FDR.
out <- "analysis"
read_tab <- function(p) read.delim(if (grepl("gz$", p)) gzfile(p) else p, check.names = FALSE)
rows <- list()
raw_p_columns <- function(columns) {
  key <- tolower(gsub(".", "_", columns, fixed=TRUE))
  columns[grepl("(^|_)(pvalue|pval|p_value|p_val|p)(_|$)|_tail$",key) &
          !grepl("^(n_|reported_|p_reproduction_)|fdr|adj|min_attainable|minimum",key)]
}
add <- function(path, family, role, correction, splitby = NULL) {
  # A compact retained release uses the canonical file.  A genuine raw rebuild
  # must inventory its new producer output rather than the archived baseline.
  if (path %in% unname(analysis_input_aliases)) {
    legacy <- names(analysis_input_aliases)[match(path, unname(analysis_input_aliases))]
    actual <- analysis_input(legacy)
    path <- substring(actual, nchar(normalizePath(".")) + 2L)
  }
  x <- read_tab(path)
  groups <- if (is.null(splitby)) list(all = x) else split(x, x[[splitby]])
  for (g in names(groups)) {
    z <- groups[[g]]
    pcols <- raw_p_columns(names(z))
    fcols <- names(z)[grepl("FDR|adj.P.Val|p_val_adj|^BH_", names(z))]
    rows[[length(rows) + 1]] <<- data.frame(family = if (is.null(splitby)) family else paste(family, g, sep = ":"),
      n_endpoints_or_retained_entries = nrow(z), n_reported_P_columns = length(pcols),
      reported_P_columns = paste(pcols, collapse = ";"), correction = correction,
      reported_FDR_columns = paste(fcols, collapse = ";"), role = role, source = path)
  }
}
add("results/tables/pseudobulk/scrna_paired_edgeR_full.tsv.gz", "discovery_lineage_gene_effects", "exploratory_n2_P_values_not_population_inference", "BH_within_each_lineage_not_used_as_discovery_evidence", "cell_type")
add("results/tables/programs/scrna_pathway_effect_summary.tsv", "discovery_pathway_effects", "descriptive_n2_no_tests", "none", "cell_type")
tf <- read_tab("results/tables/programs/scrna_gtrd_tf_target_scores_top_variable.tsv")
tf <- unique(tf[, c("cell_type", "tf_target_set")])
for (ct in unique(tf$cell_type)) rows[[length(rows) + 1]] <- data.frame(family = paste("discovery_GTRD_targets", ct, sep = ":"),
  n_endpoints_or_retained_entries = sum(tf$cell_type == ct), n_reported_P_columns = 0,
  reported_P_columns = "", correction = "none", reported_FDR_columns = "", role = "descriptive_top100_variance_selected_target_sets",
  source = "results/tables/programs/scrna_gtrd_tf_target_scores_top_variable.tsv")
comm <- read_tab("results/tables/communication/cross_method_donor_aware_edges.tsv")
for (stage in c("joined_communication_candidates", "retained_communication_edges")) {
  nr <- if (stage == "joined_communication_candidates") nrow(comm) else sum(comm$retained)
  rows[[length(rows) + 1]] <- data.frame(family = stage, n_endpoints_or_retained_entries = nr,
    n_reported_P_columns = 0, reported_P_columns = "", correction = "none", reported_FDR_columns = "",
    role = "descriptive_expression_supported_edges_no_population_tests",
    source = "results/tables/communication/cross_method_donor_aware_edges.tsv")
}
add("analysis/bulk_gate_inputs/original_A_gene_results.tsv.gz", "GSE254300_AIS_genes", "bulk_gene_level", "BH_across_tested_genes")
add("analysis/bulk_gate_inputs/original_B_gene_results.tsv.gz", "GSE311180_AIS_genes", "bulk_gene_level", "BH_across_tested_genes")
add("results/tables/meta_validation_requant/gene_cross_cohort_replication.tsv.gz", "two_cohort_Stouffer_genes", "exploratory_synthesis_not_independent_replication_screen", "BH_across_common_genes")
add("results/tables/disease_control_gse254300_cs/GSE254300_CS_paired_edgeR_full.tsv.gz", "GSE254300_CS_genes", "disease_control_gene_level", "BH_across_tested_genes")
add("results/tables/disease_control_gse254300_cs/GSE254300_AIS_by_side_interaction_edgeR_full.tsv.gz", "GSE254300_AIS_CS_gene_interactions", "diagnosis_by_side_gene_level", "BH_across_tested_genes")
add("results/tables/bulk/frozen_pathway_paired_tests.tsv", "GSE254300_AIS_frozen_pathways", "bulk_pathway_projection", "BH_across_all_retained_Hallmark_and_Reactome_pathways_t_test")
add("results/tables/bulk_gse311180_requant/frozen_pathway_paired_tests.tsv", "GSE311180_AIS_frozen_pathways", "bulk_pathway_projection", "BH_separately_across_all_pathways_for_t_and_exact_tests")
add("results/tables/disease_control_gse254300_cs/frozen_pathway_within_cohort_tests.tsv", "GSE254300_joint_frozen_pathways", "joint_normalized_within_diagnosis_pathway_projection", "BH_across_all_pathways_within_each_diagnosis_exact_tests", "cohort")
add("results/tables/disease_control_gse254300_cs/frozen_pathway_AIS_CS_interaction_tests.tsv", "GSE254300_frozen_pathway_interactions", "diagnosis_by_side_pathway_projection", "BH_across_all_pathway_interactions_exact_tests")
add("results/tables/adjacent_tissue_gse324067/GSE324067_osteoblast_paired_edgeR.tsv.gz", "GSE324067_osteoblast_paired_genes", "adjacent_tissue_context", "BH_across_tested_genes")
add("results/tables/adjacent_tissue_gse324067/GSE324067_three_site_edgeR.tsv.gz", "GSE324067_osteoblast_three_site_genes", "adjacent_tissue_context", "BH_across_tested_genes_separately_per_contrast", "contrast")
add("results/tables/adjacent_tissue_gse324067/frozen_pathway_paired_tests.tsv", "GSE324067_frozen_pathways", "adjacent_tissue_context", "BH_separately_across_all_pathways_for_t_and_exact_tests")
add("results/tables/reference_gse17503/GSE17503_paired_limma_full.tsv.gz", "GSE17503_biopsy_culture_genes", "culture_reference_context_not_side_validation", "BH_across_tested_genes")
add("results/tables/reference_gse17503/cell_type_marker_paired_tests.tsv", "GSE17503_lineage_marker_scores", "culture_reference_context_not_side_validation", "BH_separately_across_marker_scores_for_t_and_exact_tests")
add("results/tables/robustness/pseudobulk_vs_cell_mixed_model.tsv", "discovery_cell_mixed_model_diagnostics", "effect_sign_diagnostic_singular_models_not_population_inference", "no_population_inferential_use")
add("analysis/bulk_patient_effect_statistics.tsv", "revision_original_score_cohorts_and_pool", "same_frozen_endpoint_alternative_tests", "nominal_no_overall_FDR_across_cohorts_or_test_types")
add("analysis/bulk_pooled_leave_one_patient_out.tsv", "revision_pooled_patient_omissions", "sensitivity_same_endpoint_not_independent_validation", "nominal_no_overall_FDR")
add("analysis/bulk_age20_and_site_exclusion.tsv", "revision_age_site_score_sensitivity", "locked_or_restandardized_frozen_endpoint_sensitivity", "nominal_no_overall_FDR")
add("analysis/bulk_cohort_stratified_estimates.tsv", "revision_cohort_factor_models", "two_explicit_estimands", "nominal_no_overall_FDR")
add("analysis/bulk_meta_analysis_heterogeneity.tsv", "revision_meta_models", "two_cohorts_Hartung_Knapp_and_fixed_effect_sensitivity", "nominal_no_overall_FDR")
add("analysis/bulk_interaction_uncertainty.tsv", "revision_score_interaction", "diagnosis_by_side_patient_effect", "nominal_no_overall_FDR")
add("analysis/bulk_gene_concordance_patient_permutation_summary.tsv", "revision_genomewide_concordance", "two_genomewide_patient_side_statistics", "nominal_no_overall_FDR")
if(file.exists("analysis/bulk_gene_concordance_edgeR_permutation_summary.tsv"))
  add("analysis/bulk_gene_concordance_edgeR_permutation_summary.tsv", "revision_original_edgeR_genomewide_concordance", "two_original_estimator_genomewide_patient_side_statistics", "nominal_no_overall_FDR")
add("analysis/bulk_RNA_pairing_gate_exclusion.tsv", "revision_RNA_pairing_gate_exclusion", "post_hoc_locked_or_restandardized_fixed_endpoint_sensitivity", "nominal_no_overall_FDR")
add("analysis/bulk_followup_ADIPOQ_proxy_models.tsv", "followup_ADIPOQ_expression_proxy_intercepts", "six_post_hoc_target_by_cohort_intercepts_not_measured_adiposity", "BH_across_six_intercepts_no_study_wide_FDR")
add("analysis/bulk_followup_available_clinical_modifiers.tsv", "followup_available_age_Cobb_modifiers", "four_post_hoc_slopes_plus_two_centered_mean_intercepts", "BH_across_four_slopes_centered_intercepts_nominal_no_study_wide_FDR")
# Additional retained P-value artifacts omitted by the earlier partial inventory.
# These blocks overlap: selected subsets, alternative estimators and archived
# recipes are explicitly identified and must not be summed as independent tests.
add("results/tables/annotation/cluster_differential_markers.tsv.gz", "annotation_cluster_markers", "cell_level_annotation_only_not_patient_side_inference", "Seurat_Bonferroni_marker_adjustment_not_condition_FDR", "cluster")
add("results/tables/annotation/top15_cluster_markers.tsv", "annotation_selected_top15_markers", "selected_display_subset_of_annotation_cluster_markers_not_new_tests", "same_adjusted_values_as_full_annotation_marker_table")
add("results/tables/bulk/multicellular_axis_model_sensitivity.tsv", "original_first_cohort_score_models", "overlapping_fits_of_same_score_not_independent_endpoints", "nominal_no_overall_FDR")
add("results/tables/bulk_gse311180_requant/frozen_axis_test.tsv", "original_second_cohort_score", "duplicate_original_score_also_in_revision_statistics", "nominal_no_overall_FDR")
add("results/tables/bulk_gse311180_requant/frozen_axis_leave_one_patient_out.tsv", "original_second_cohort_patient_omissions", "original_within_cohort_sensitivity", "nominal_no_overall_FDR")
add("results/tables/bulk_gse311180_requant/frozen_axis_remove_top1pct.tsv", "original_second_cohort_high_expression_exclusion", "post_hoc_score_sensitivity", "nominal_no_overall_FDR")
add("results/tables/bulk_gse311180_requant/frozen_axis_random_signature_test.tsv", "original_second_cohort_matched_signature", "archived_recipe_also_in_revision_recipe_audit", "nominal_empirical_tail_no_overall_FDR")
add("results/tables/disease_control_gse254300_cs/frozen_axis_within_cohort_tests.tsv", "original_joint_score_by_diagnosis", "joint_scaling_distinct_from_AIS_only_projection", "nominal_no_overall_FDR")
add("results/tables/disease_control_gse254300_cs/frozen_axis_model_sensitivity.tsv", "original_joint_score_models", "within_diagnosis_covariate_sensitivities", "nominal_no_overall_FDR", "cohort")
add("results/tables/disease_control_gse254300_cs/frozen_axis_AIS_CS_interaction_test.tsv", "original_score_diagnosis_interaction", "duplicate_endpoint_also_in_revision_interaction_statistics", "nominal_no_overall_FDR")
add("results/tables/meta_validation_requant/frozen_axis_12_patient_test.tsv", "original_pooled_score", "duplicate_endpoint_also_in_revision_statistics", "nominal_no_overall_FDR")
add("results/tables/robustness/bulk_axis_exact_side_label_permutation.tsv", "original_first_cohort_exact_score", "duplicate_endpoint_also_in_revision_statistics", "nominal_no_overall_FDR")
add("results/tables/robustness/bulk_axis_matched_random_gene_sets_summary.tsv", "original_first_cohort_matched_signature", "archived_124_slot_recipe_not_corrected_122_slot_recipe", "nominal_empirical_tail_no_overall_FDR")
add("results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_test.tsv", "original_osteoblast_score", "adjacent_tissue_not_bilateral_muscle_validation", "nominal_no_overall_FDR")
add("results/tables/adjacent_tissue_gse324067/frozen_muscle_axis_leave_one_patient_out.tsv", "original_osteoblast_patient_omissions", "adjacent_tissue_score_sensitivity", "nominal_no_overall_FDR")
add("results/tables/published_summary_gse311180/frozen_axis_random_set_test.tsv", "published_summary_signature_overlap", "published_only_gene_list_overlap_not_patient_expression_validation", "nominal_empirical_tail_no_overall_FDR")
add("results/tables/published_summary_gse311180/frozen_pathway_random_set_tests.tsv", "published_summary_pathway_overlaps", "published_only_gene_list_overlap_not_pathway_activation", "nominal_empirical_tails_no_overall_FDR")
add("analysis/bulk_cohort_mean_difference.tsv", "revision_between_cohort_score_difference", "processing_heterogeneity_diagnostic", "nominal_no_overall_FDR")
add("analysis/bulk_cohort_random_null_recipe_audit.tsv", "revision_separate_matched_signature_recipes", "archived_and_corrected_recipes_overlapping_benchmarks", "nominal_empirical_tails_no_overall_FDR")
add("analysis/bulk_pooled_matched_random_signature_summary.tsv", "revision_joint_matched_signature_recipes", "same_random_gene_sets_projected_in_two_cohorts_and_pool", "nominal_empirical_tails_no_overall_FDR")
add("analysis/bulk_gene_heterogeneity_logCPM_descriptive.tsv", "revision_gene_heterogeneity", "two_cohort_descriptive_mean_logCPM_estimator_not_original_QL_meta_analysis", "BH_across_common_genes")
add("analysis/bulk_interaction_within_diagnosis_statistics.tsv", "revision_joint_score_within_diagnosis", "alternative_tests_of_original_joint_scaled_scores", "nominal_no_overall_FDR")
add("analysis/signature_bulk_sensitivity_tests.tsv", "revision_signature_subset_scores", "five_lineage_only_tests_BH_other_subsets_nominal_out_of_gate_subsets_descriptive", "BH_only_across_five_lineage_subscores_within_each_cohort_or_pool_other_rows_nominal", "cohort")
add("analysis/signature_literature_nominated_bulk_effects.tsv", "revision_literature_gene_lookup", "selected_lookup_of_original_gene_tests_not_new_seven_gene_family", "original_genomewide_gene_family_FDR_no_recorrection")
add("analysis/signature_pipeline_flip_projection_tests.tsv", "revision_discovery_flip_bulk_projection_tests", "four_reselected_signatures_projected_in_each_cohort_and_pool", "nominal_no_overall_FDR", "cohort")
add("analysis/signature_pipeline_flip_null_summary.tsv", "revision_discovery_flip_pipeline_null", "four_configurations_two_unsigned_selections_no_gene_or_cell_independence", "nominal_exact_pipeline_tail_no_overall_FDR")
# Bounded biological revision: these remain retrospective and overlap across normalisation.
add("analysis/bulk_gate_genomewide_exact_summary.tsv", "RNA_gate_fixed_genomewide_endpoints", "post_hoc_fixed_universe_patient_side_null", "BH_two_endpoints_within_normalization_and_full_four_endpoint_sensitivity_family", "normalization")
for(norm in c("locked_original_TMM","retained_library_reTMM"))
  add(paste0("analysis/bulk_gate_observed_genes_",norm,".tsv.gz"), paste0("RNA_gate_full_genes_",norm), "fixed_17961_gene_refit_no_reselection", "BH_across_original_17961_gene_universe")
add("analysis/bulk_gate_ANGPT1_LEP_fixed_candidate_results.tsv", "RNA_gate_prior_fixed_candidates", "retrospective_two_gene_lookup_full_gene_FDR_also_reported", "original_full_gene_FDR_plus_explicit_BH_two_candidates_and_four_candidate_normalization_sensitivity_family", "normalization")
add("analysis/bulk_gate_candidate_patient_exact.tsv", "RNA_gate_candidate_patient_resolution", "later_added_retrospective_resolution_check_same_two_candidates", "BH_full_four_candidate_normalization_family_no_study_wide_correction", "normalization")
add("analysis/biological_Cheng2026_fixed_target_public_summary.tsv", "independent_published_target_availability", "no_available_fitted_estimate_or_test_P_all_NA_not_measured_null", "no_test_calculated_not_independent_validation")
add("results/tables/communication/LIANA_NATMI_CellPhoneDB_all.tsv.gz", "communication_method_score_P_values", "algorithmic_cellphone_method_scores_not_donor_side_inference", "native_method_scores_no_patient_level_multiple_testing_or_independent_evidence", "method")
add("analysis/communication_native_precision_boundaries.tsv", "communication_native_P_precision_source",
    "precision_source_copy_of_37_existing_CellPhoneDB_method_P_values_not_new_tests_or_family",
    "inherits_native_method_support_rule_no_new_inference_or_multiple_testing_family")
# Joined coefficient-comparison tables repeat constituent tests and original FDR.
add("results/tables/adjacent_tissue_gse324067/muscle_osteoblast_effect_comparison.tsv.gz", "cross_tissue_gene_lookup", "joined_existing_gene_tests_descriptive_cross_tissue_comparison_not_new_family", "inherits_separate_original_tissue_gene_family_FDR_no_combined_correction")
add("results/tables/disease_control_gse254300_cs/GSE254300_AIS_CS_gene_effect_comparison.tsv.gz", "cross_diagnosis_gene_lookup", "joined_existing_AIS_CS_interaction_tests_not_new_family", "inherits_each_original_gene_family_FDR_no_lookup_recorrection")
add("results/tables/published_summary_gse311180/GSE254300_GSE311180_effect_comparison.tsv.gz", "published_summary_gene_lookup", "joined_original_AIS_tests_and_published_summary_not_patient_expression_validation", "inherits_original_AIS_gene_family_and_published_source_FDR_no_new_combined_family")
# Figure-source exports repeat already reported tests. Include their artifacts
# for traceability, without treating copies or old displays as new test families.
for (directory in c("results/figures/source_data", "analysis/tables")) {
  for (path in sort(list.files(directory, pattern="\\.tsv(\\.gz)?$", full.names=TRUE))) {
    if (length(raw_p_columns(names(read_tab(path)))))
      add(path, paste0("display_source:",basename(path)),
          "duplicate_or_archived_display_of_existing_tests_not_independent_family",
          "inherits_original_test_correction_no_new_display_family_or_study_wide_FDR")
  }
}
inventory <- do.call(rbind, rows)
write.table(inventory, file.path(out, "bulk_test_family_inventory.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(inventory), "retained original/revision analysis blocks inventoried.\n")
