# Scoliosis muscle transcriptomic reanalysis

The project analyses public bilateral paraspinal-muscle expression data from
PRJNA722100, GSE254300 and GSE311180, with GSE324067 and GSE17503 as
role-restricted tissue and culture controls. Biological sample size is the
number of patients, not the number of cells or genes.

## Files

- `code/`: scientific analysis and figure scripts, configuration and environments.
- `data/processed_data_*.zip`: processed counts, sample metadata, annotation evidence,
  pseudobulk results, published numerical results and figure-source tables as
  ordinary TSV files. All logical table paths and their archive members are indexed by `data/CONTENTS.tsv`; exact duplicate tables share a single member.
- `data/results/tables/bulk/frozen_multicellular_signature.tsv`: selected
  lineage–gene entries and discovery weights used by scoring scripts.
- `MANIFEST.tsv`: file sizes and SHA-256 checksums for the release.

The data archives group related tables for distribution without changing their
scientific content. `data/CONTENTS.tsv` permits each table to be identified and
checked independently. No external reference-output cache is required.

## Run

Install the dependencies in `code/environment/environment.yml`, then run from
the project directory:

```sh
python code/run.py
```

This checks the source files, extracts scientific tables to `results/analysis/`,
recomputes effective-length-aware bulk models, shared-gene projections and
joint/conditional patient-side permutations, and rebuilds signature-level cell
expression summaries and the marker figure. Calculated results are written to
`results/analysis/computed/`

