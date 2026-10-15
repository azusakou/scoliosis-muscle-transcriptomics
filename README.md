# Scoliosis muscle transcriptomic reanalysis

- `code/`: analysis programs, tests, configuration and recorded environments.
- `data/`: processed inputs, source metadata and frozen expected results.

With the recorded environment installed, run from the extracted repository root:

```sh
python code/run.py
```

The pipeline validates input data, reconstructs cell-level datasets, refits statistical models, and regenerates figures. Dependencies and recorded environments are provided in code/environment/.

Raw sequencing data are available through GEO/SRA (PRJNA722100, GSE254300, GSE311180, GSE324067, GSE17503). 




