# Reproducing the analysis

This repository contains scripts, not input datasets. Download inputs from the original sources listed in `DATA_SOURCES.md`. Do not place downloaded data or generated result files under version control.

## Python workflows

Run Python scripts from GBM_hypoxia_glutamine-GitHub/, the directory containing python/ and R/. Scripts use repository-relative `data/` and `results/` defaults. The expected input tree includes:

```text
data/
  cohort_audit/       # source evidence required by the manifest builder
  lin/                # Lin cohort .h5ad input
  greenwald_visium/   # extracted Greenwald Visium sections
  ravi_visium/        # extracted Ravi Visium sections
  gse243848/          # normalized GSC expression table
  gse300765/          # GEO series matrix
results/
```

Main stages are organized under `python/00_manifest` through `python/05_figures`: cohort manifest preparation; Lin and spatial analyses; GSC and U87MG analyses; supplementary/statistical recomputations; and figure-data/figure generation. Some scripts consume intermediate files from earlier stages. Check each script's input/output constants before running it, and run stages in order. The source inputs and intermediate tables are not included in this repository.

Optional environment variables can redirect paths:

- `GBM_DATA_DIR`: root for downloaded inputs (default `data/`).
- `GBM_OUTPUT_DIR`: root for generated results (default `results/`).
- `GBM_LIN_H5AD`, `GBM_GREENWALD_DIR`, `GBM_RAVI_DIR`, `GBM_GSC_FILE`, `GBM_U87_FILE`: explicit input paths.
- `GBM_COHORT_MANIFEST`: explicit cohort manifest path.
- `GBM_U87_DIR`, `GBM_U87_SCORES`: U87 annotation/output and score paths used by the relevant scripts.

Install Python 3.11 and the dependencies required by the scripts you intend to run. The key-version file is a record of selected versions, not a complete install specification.

## R workflows

Run from the `R/` directory so relative paths resolve as expected. Put the required inputs under `R/data/`. `00_setup.R` defines shared settings and helper functions. Run `01_gateA_lin.R` through `04_gateC_u87.R` before `05_figures_main.R` and `06_figures_supp.R`.

The workflows use R 4.4.3 and packages including data.table, GEOquery, ggplot2, ggrepel, limma, Matrix, matrixStats, patchwork, RColorBrewer, scales, Seurat/SeuratObject 5, SingleCellExperiment, zellkonverter, hdf5r, svglite, and reshape2. Some packages require system libraries, including HDF5, zlib, libuv, and FreeType. Install package versions compatible with your operating system and R version.

## Outputs

Scripts write analysis tables and figures beneath their configured output directories. Outputs are reproducible products of the input data and scripts and are not committed here. Consult the script headers for precise filenames and any stage-specific requirements.
