# HLA-Splitter Analysis

Code accompanying **Label-free HLA-guided demultiplexing enables single-cell analysis of scarce clinical samples to reveal disease-specific immune signatures**.

## Structure

| Directory | Analysis | Supplementary Figures |
| --- | --- | --- |
| `Figure1_benchmark` | Synthetic pools, accuracy/recovery/runtime, HLA typing and mixture validation | S1-S3 |
| `Figure2_CSF` | Immune composition, myeloid programs and T-cell analysis | S4-S5 |
| `Figure3_BALF` | T-cell states, modules, signatures and TCR expansion | S6 |
| `Figure4_urine_immune` | Immune composition, rejection modules and cross-cohort signatures | S7-S8 |
| `Figure5_urine_Tcell_CXCR3` | CXCR3, TCR, ligand-receptor analysis and GSE189536 validation | S9-S10 |
| `00_Shared` | Reusable plotting, enrichment and path helpers | Shared |

See directory READMEs for individual analyses. Numbers organize the files, not an automatic execution order. Upstream scripts require additional raw inputs; reuse annotated objects and source tables for figure regeneration.

Manual BAM-based benchmarks are in [Figure1_benchmark/BAM_demultiplex](Figure1_benchmark/BAM_demultiplex/README.md). They are separate from the R figure workflows.

## Environment

Reference R version: **4.5.3**. Install packages needed by the selected analysis:

```r
install.packages(c("Seurat", "harmony", "Matrix", "dplyr", "tidyr", "tibble",
  "readr", "readxl", "stringr", "purrr", "data.table", "ggplot2", "ggpubr",
  "patchwork", "cowplot", "ggrepel", "ggsci", "scales", "pheatmap", "RColorBrewer",
  "circlize", "ragg", "pROC", "RANN", "Rtsne", "enrichR", "BiocManager"))
BiocManager::install(c("AUCell", "ComplexHeatmap", "clusterProfiler", "GSVA",
  "limma", "MAST", "scRepertoire"), ask = FALSE, update = FALSE)
```

Additional dependencies appear in each script. Optional CellChat workflows need their respective packages. Exports use Arial and Cairo/ragg where specified.

## Data and Paths

Keep the processed-data bundle beside `Code`:

```text
Revison/
|-- Code/
`-- Data/
    |-- synthetic/
    |-- pbmc30/
    |-- pbmc4/
    |-- hash10x/
    |-- csf/
    |-- balf/
    |-- urine/
    |-- gse189536/
    |-- references/
    `-- data_manifest.csv
```

Dataset folders contain annotated objects, metadata and required source tables. `data_manifest.csv` records input checksums. Raw reads and general intermediate files are excluded. Paths are written explicitly in the individual R scripts; edit each input, output and helper-script path to match your files. The supplied paths point to the current local project layout.

## Manual Use

Open the required R script in R/RStudio, review its paths and run the relevant sections manually. Load any required saved objects and complete annotation checkpoints before proceeding. Some upstream scripts require an initialized analysis session rather than a clean session.

There is no shared execution entry point, automatic path configuration or automatic whole-figure workflow. Existing statistical settings and named palettes remain in the individual scripts.

For Linux BAM-based benchmarking, edit `Figure1_benchmark/BAM_demultiplex/config.sh` and invoke the required stage scripts directly as described in its README.
