# Retained data definitions and roles

All expression comparisons use the patient/donor as the independent unit.
Cells and genes do not increase the biological sample size. An omitted cell
lineage with a verified zero count is represented explicitly as zero in every
library's composition grid; unavailable expression features are not imputed.

| Field or family | Meaning |
| --- | --- |
| `sample_id`, `processed_sample_id`, `pb_id` | Deposit/processing identifiers, not proof of unique specimen identity |
| `patient_id`, `pair_id` | Publicly annotated pairing; `pair_id` includes diagnosis where needed |
| `side` | Concave/convex deposit label; this does not independently verify anatomical left/right |
| `paired_difference` | Concave score minus convex score; biopsy-minus-culture in GSE17503 tables |
| `weight` | Discovery-donor direction, averaged across lineage entries for a repeated gene; zero conflicts excluded |
| Projection `score` | Sum(weight × within-cohort gene z-score) / sum(abs(weight)); sample SD used; flat rows become zero |
| `logFC` | Patient-blocked edgeR or stated limma coefficient, not the average paired TMM-logCPM statistic |
| `paired_logCPM_difference` | Donor concave-minus-convex TMM log2-CPM contrast |
| `PValue`, `p_value`, `t_test_p` | Raw P value under the named test; exploratory where specified |
| `FDR`, `FDR_exact`, other FDR columns | Benjamini–Hochberg adjustment within the explicitly named family; partial overlap tables retain original full-family q values |
| `exact_sign_flip_p` | Complete patient-side sign-flip test; finite P resolution depends on the number of pairs |
| Empirical random-signature P | (1 + number of absolute null effects at least as large as observed) / (1 + draws) |
| `Freq`, `n_cells`, `final_singlets` | Retained annotated-cell counts; composition fractions sum within library |
| `proportion` | Fraction 0–1; NNLS/ridge-NNLS weights are descriptive estimates and lack external accuracy validation |
| `percent_expressing` | Despite the historical name, fraction 0–1, not a percentage 0–100 |
| `mean_normalized_expression` | AverageExpression's mean after expm1 of the data layer; formerly mislabeled `mean_log_expression` |
| `average_log_expression` | Direct average of the log-normalized layer in marker specificity tables; distinct from AverageExpression |
| `side_semipartial_r2` | Sequential side sum of squares after patient adjustment divided by total component sum of squares; formerly mislabeled partial R² |
| `sample_r2` | R² of the component model with sample ID as factor |
| `inferential_status`, `interpretation_scope` | Required limits on biological and statistical interpretation |
| `NA` / `not retained/detectable` | Missing/unavailable value, not a numerical zero |

`results/tables` contains original retained analysis outputs. The archived
1,000-draw GSE254300 random null used 124 slots, including two zero-weight
conflicts; inference in the revised paper uses the corrected 122-slot
10,000-draw null in `analysis`. The archived null is preserved to expose the
change rather than replace history. GSE254300 high-expression sensitivity
counts now describe nonzero scoring genes: 122 before trimming, 113 after;
its score values did not change. The GSE311180 trimming recipe is separate.

`analysis/*.tsv` and `analysis/tables` contain completed reviewer
statistics and revised figure source data. `results/figures/source_data`
contains historical panel source tables; escaped multiline labels are valid
TSV fields. `metadata/audit` contains source maps,
version/checksum records and role restrictions, including historical local
paths for provenance. Historical paths are not executable defaults.

GSE324067 osteoblasts, GSE17503 biopsy/culture and published-only summaries
remain role-restricted context; they are not pooled into bilateral-muscle
validation. Subset tables' retained adjusted P values must not be read as
adjusted over only the subset rows.

Real-cell expression audit: `assay=RNA` means uncorrected counts;
`assay=SoupX` means corrected counts in the same retained cells/labels.
`detection_fraction` is 0–1; `n_detected` counts cells with UMI>0.
`mean_UMI_per_cell` includes zero-count cells, while
`mean_UMI_detected_cells` is undefined when no cell is positive.
`mean_cp10k` averages UMI/(all-gene per-cell total)×10000.
`*_log1p_cp10k` uses the natural logarithm after this normalization; the 90th
percentile uses R's type-7 quantile. `n_cells=0` groups have undefined expression
and detection, not estimated zero expression. The full summary covers158genes;
the 296,920-row source table and independent validator cover ten nominated genes
in all14,846cells and both assays. FigureS7 depicts these ten genes. Source,
script and software hashes are in the completed real-cell manifest.


RNA-failure diagnostics (all post hoc; primary verdicts unchanged):

- `audit_RNA_pairing_loss_decomposition.tsv`: seven expected pairs; integer
  het/hom and opposite-homozygote counts, full-precision relatedness, additive
  contributions to its deficit, and direction of het/hom discordance.
- `audit_RNA_fingerprint_site_coverage.tsv`: 14 libraries and the common 16,469
  autosomal-marker denominator; zero/covered/callable counts, rejection reasons,
  depth medians and fraction of depth in the 165 most covered markers.
- `audit_RNA_pairing_depth_strata.tsv` / `audit_RNA_pairing_chromosome_strata.tsv`: shared-call
  partitions defined by minimum per-library depth or 22 autosomes. Additive
  integer fields sum to the primary pair; relatedness values do not sum.
- `audit_RNA_pairing_depth_AB_grid.tsv`: 12 settings × 91 comparisons; altered
  calling settings are not altered primary interpretation thresholds.
- `audit_RNA_pairing_directional_dropout.tsv`: 14 directional summaries;
  homozygous-side depth bins, zero-minor-allele counts and retained reference/
  alternate allele counts. Alphabetical sketch count order is explicitly
  mapped to VCF reference order; per-marker variants/counts are not published.
- `audit_RNA_pairing_unique_high_quality_reads.tsv`: 91 strict-pileup
  comparisons; MAPQ60/baseQ20, uncapped depth and BAQ disabled. The canonical
  table uses samtools paired-overlap adjustment; `_overlap_retained.tsv` instead
  counts both overlapping mates. Both are independent pileup conventions, not
  exact native extraction. Corresponding transition tables preserve both modes.
- `audit_RNA_read_filter_call_transitions.tsv`: 7 pairs × 4 original call
  categories × 4 strict call categories, including uncalled; each pair totals
  all 16,469 markers, preserving the same-marker denominator.
- `audit_RNA_pairing_full_read_cap50.tsv` and `audit_RNA_cap50_source_QC.tsv`:
  15 comparisons among 6 completely re-aligned libraries, including both
  nonpassing pairs and one passing control; full read counts/source/tool-bound
  checks and artifact hashes. Separate from the uniform 14-library cap20 primary.
- `audit_RNA_cap50_call_transitions.tsv`: 3 pairs × 4 primary categories ×
  4 cap50 categories, 48 records; each pair totals all 16,469 markers.
- `audit_RNA_CHARR_proxy.tsv` / `audit_RNA_paired_MLE_proxy.tsv`: 28 sample
  records / 364 ordered comparisons at depth 15/30. Raw native method outputs,
  explicitly uncalibrated for RNA; values must not be interpreted as admixture
  percentages, DNA identity, anatomical labels or a diagnostic exclusion test.

- `audit_RNA_discordant_fragment_coordinate_diversity.tsv`: all14 directional
  summaries over the1,520 original het/hom-discordant loci across seven pairs.
  Primary proper pairs, MAPQ60 and individual baseQ20; mate overlap is collapsed
  by QNAME/template-coordinate entries and ambiguous/other allele entries are
  excluded. Coordinates use chromosome,min(POS,PNEXT),abs(TLEN). Medians include
  zero-coverage loci, and counts distinguish original loci from qualifying
  proper templates. These are alignment-diversity proxies, not UMI molecule
  counts or calibrated PCR duplicate rates. No read names/coordinates are released.
- `audit_RNA_cap50_alignment_attempts.tsv`: source-bound terminal receipts for
  all six requested libraries, including failed stages; an unsuccessful native
  attempt is not a failed RNA identity test or a verified pair comparison.


Compact cell replay inputs under `data/retained/` retain the nonzero integer
counts for the fixed union of the 158-gene expression audit and 18 biological
targets (165 unique genes). `cell_assay_features.tsv.gz` lists the complete
assay feature universe, including features with zero target counts.
`cell_target_availability.tsv` distinguishes absent features from available
zero-count genes. Normalization always uses the original all-gene cell totals
in the annotated-cell metadata. Fresh replay reports bind the sparse inputs;
they do not reuse historical RDS object hashes as their own input identity.

`analysis/communication_native_precision_boundaries.tsv` contains 37 native
CPDB P values whose printed TSV spelling is 0.05 while the source double is
0.050000000000000044. It restores precision for the original support rule,
not a new threshold or family. The source manifest records the actual cache
hashes and complete-table field matching. Collapsed complex keys can have
multiple legacy communication rows; validation preserves this multiset and
does not count repeated rows as independent biological interactions.

`data/retained/mixed_model_extra/` adds genuine SoupX SLN counts for all
183 retained SMC-like cells, including 172 explicit zero counts. The feature
availability and reconstructed-object/exporter hashes are recorded separately
from the unchanged 165-gene export. This closes the one missing gene required
by the original 50 mixed models; no unexported feature is silently set to zero.

`data/retained/gse254300_specificity/` contains 14,967 fixed joint-universe
genes for the ten CS libraries, and only the 671 joint-universe genes missing
from the existing 14,296-gene AIS count matrix. Combining those genuine AIS
extras with the existing gate counts recreates the joint matrix without
redistributing the same AIS counts. The source export verifies the full
original count XLSX and original filter universes before writing this subset.
The refitter uses the original 14,450-gene CS and 14,967-gene interaction
models; it does not select a new favorable testing universe.
