# Data sources and input requirements

This code-only candidate does not redistribute raw data, cohort manifests, patient-level derived tables, or figure source-data CSVs. Obtain inputs from the original repositories and follow their record-specific terms.

| Cohort | Source | Role in analysis |
|---|---|---|
| Lin scRNA-seq | Zenodo record 17117252; raw human sequence data under NGDC HRA004942 | Discovery cohort |
| Greenwald Visium | GEO GSE237183; processed dataset DOI 10.5281/zenodo.8105466 | Spatial cohort |
| Ravi Visium | Zenodo record 6972429 / Dryad DOI 10.5061/dryad.h70rxwdmj | Independent spatial cohort |
| GSC hypoxia | GEO GSE243848 | 24-hour hypoxia perturbation in patient-derived cell lines |
| U87MG | GEO GSE300765 | Established cell-line boundary condition |

The following GEO input files are used by the corresponding scripts:

| Input | Official URL | SHA256 |
|---|---|---|
| GSE243848 normalized counts | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243848/suppl/GSE243848_Original_Counts.deseq2_normalized_GSC.txt.gz | `8a57865ae40f49e80c3ec24db86f2eb21031c78b406e102a6d5e96effb6e22b8` |
| GSE300765 series matrix | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE300nnn/GSE300765/matrix/GSE300765_series_matrix.txt.gz | `a4c30e0c8b6018bf097d7293cb6bade47164cf8da8a0d3a3039e97258f0e88a5` |

Those files are not included here; download them directly from GEO before running the corresponding workflows. Other raw matrix inputs and locally generated intermediate tables must also be obtained or regenerated as described by each script. The original GSE300765 metadata contains inconsistent platform-version wording; the verified mapping used the recorded probe IDs and GPL10558 annotation.

The article and this file cite the data sources; citations do not replace source-specific license or data-use terms. This repository's eventual code license will apply to project-authored code only.

