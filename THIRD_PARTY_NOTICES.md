# License scope and external resources

The MIT license applies to original analysis code and its original documentation.
It does not relicense public-deposit data, publication supplements, gene-set
resources, third-party R/Python packages, Seurat references or GENCODE annotations.
No raw sequencing reads, annotated cell objects, source article full texts or
licensed GMT files are redistributed in this release. Retained result tables
are provided as scientific source data with deposit/publication provenance;
this notice does not assign a new license to underlying third-party material.

Use the original providers' terms for PRJNA722100, GSE254300,
GSE311180/PRJNA1368770, GSE324067, GSE17503, GENCODE v49, MSigDB, GTRD,
celldex and LIANA/Consensus resources. Source/version/checksum records are in
`data/metadata/audit/`; package versions are in
`code/environment/environment.lock.tsv` and the recorded runtime table in `data/expected_results.json.gz`.
The reference-version manifest is a record of the original runtime, not a
Conda/renv solver lockfile. External tools retain their respective licenses.

LIANA was installed from saezlab/liana at the recorded full commit
6cab46c54234f861ea176c3de77c4b8aa45ecb3d, version 0.1.14. It is fetched by the
explicit installation script; it is not vendored into this release.

The code license text is the standard [SPDX MIT text](https://spdx.org/licenses/MIT.html).
Repository publication, release tagging and an optional archival DOI remain
with the author. No existing repository URL or archival identifier is asserted.

The separate numerical diagnostic `code/analysis/scripts/probe_TMM_rank_ties.R` adapts edgeR 4.4.2 TMM code and retains edgeR's GPL (>= 2) terms (SPDX: GPL-2.0-or-later), rather than the MIT license for original analysis code. The GPL version 2 text is supplied in `code/licenses/GPL-2.0-or-later.txt`. Its comments identify the upstream package/function; installed edgeR retains its original authorship and package notices. This diagnostic license does not relicense underlying research data.
