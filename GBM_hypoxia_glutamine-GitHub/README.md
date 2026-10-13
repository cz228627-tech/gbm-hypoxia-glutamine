# Hypoxia–glutamine analysis code in glioblastoma

This repository contains the Python analysis scripts and independent R statistical workflows for the study *Hypoxia-associated glutamine transcription in glioblastoma is reproducible across patients but reorganized by biological context*.

The repository is code-only. It does not contain downloaded expression data, patient-level data, cohort manifests, or generated figure/source-data tables. Obtain the inputs from the accessions listed in [DATA_SOURCES.md](DATA_SOURCES.md), and follow each source's terms.

## Contents

- `python/`: analysis scripts, arranged by cohort and workflow stage.
- `R/`: independent R analyses and figure-generation scripts.
- `RUN_GUIDE.md`: input layout, execution order, and environment variables.
- `DATA_SOURCES.md`: public data accessions and input notes.
- `requirements-python-key-versions.txt`: recorded key Python package versions; this is not a complete lockfile.
- `LICENSE`: MIT license for the project code in this repository.

## Software

The recorded analysis environment used Python 3.11.13 and the package versions listed in `requirements-python-key-versions.txt`. The R workflows were recorded with R 4.4.3, data.table 1.18.6.1, and limma 3.58.1; other R and system dependencies are described in `RUN_GUIDE.md`. These records document the analysis environment but do not fully pin every dependency.

## Reuse and citation

Please cite the associated manuscript and the original data sources. The `CITATION.cff` file provides software citation metadata. Cite the version-specific Zenodo DOI after the first public release is archived.

The MIT license covers the code authored for this project only. It does not grant rights to datasets, third-party software, or other works.

