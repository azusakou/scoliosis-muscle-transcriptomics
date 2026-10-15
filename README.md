# Scoliosis muscle transcriptomic reanalysis

- `code/`: analysis programs, tests, configuration and recorded environments.
- `data/`: processed inputs, source metadata and frozen expected results.

With the recorded environment installed, run from the extracted repository root:

```sh
python3 code/run.py
```

Outputs and check reports appear in `results/run-*/`; the 16 main/supplementary
figure PDFs are under that run's `workspace/results/figures/` and
`workspace/analysis/figures/`. A successful run ends with `PASS.json`.

The command verifies `MANIFEST.tsv`, reconstructs four cell-level expansions
from authentic sparse counts, independently joins them back to counts and full
library totals, refits the statistical models, regenerates the figures, compares
against frozen results, and runs regression/independent checks. Frozen tables
are consolidated in `data/expected_results.json.gz`; every original table cell
is retained as a string. Generated outputs never replace reference inputs.

The recorded runtime uses R 4.4.2, edgeR 4.4.2, lme4 2.0.1 and metafor 4.8.0.
Environment specifications, version records and platform locks are in
`code/environment/`. These include dependencies for the retained upstream
scientific recipes. Raw reads remain at GEO/SRA (PRJNA722100, GSE254300,
GSE311180, GSE324067 and GSE17503); this command reproduces the paper from
processed inputs, without repeating read alignment or upstream cell inference.
Data definitions are in `code/docs/DATA_DICTIONARY.md`.

The original score remains unconfirmed after excluding the two failed RNA
pairs. RNA consistency cannot authenticate anatomical side. Cell mixed models
remain diagnostics from two patients (50 models, 49 sign agreements, 38 singular
fits). Missing clinical identity and biological experiments remain limitations.
See `LICENSE` and `THIRD_PARTY_NOTICES.md` for code/data terms. Publication and
the actual repository URL/archival identifier remain with the author.
