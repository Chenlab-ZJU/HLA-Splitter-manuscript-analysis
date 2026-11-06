# HLA-Splitter Analysis

Code accompanying **Label-free HLA-guided demultiplexing enables single-cell analysis of scarce clinical samples to reveal disease-specific immune signatures**.

## Code Structure

| Directory | Main Figure | Supplementary Figures | Analysis |
| --- | --- | --- | --- |
| [Figure1_benchmark](Figure1_benchmark/README.md) | Fig. 1 | S1-S3 | Demultiplexing benchmarks and HLA/PBMC validation |
| [Figure2_CSF](Figure2_CSF/README.md) | Fig. 2 | S4-S5 | CSF immune composition, myeloid programs and T cells |
| [Figure3_BALF](Figure3_BALF/README.md) | Fig. 3 | S6 | BALF T-cell states, signatures and TCR expansion |
| [Figure4_urine_immune](Figure4_urine_immune/README.md) | Fig. 4 | S7-S8 | Urine immune composition, rejection modules and signatures |
| [Figure5_urine_Tcell_CXCR3](Figure5_urine_Tcell_CXCR3/README.md) | Fig. 5 | S9-S10 | CXCR3, TCR, communication and biopsy comparisons |
| `00_Shared` | Shared | Shared | Plotting/enrichment helpers and metadata export |

## Environment

Reference R version: **4.5.3**. Common dependencies:

```r
install.packages(c("Seurat", "harmony", "Matrix", "dplyr", "tidyr", "tibble",
  "readr", "readxl", "stringr", "purrr", "data.table", "ggplot2", "ggpubr",
  "patchwork", "cowplot", "ggrepel", "ggsci", "scales", "pheatmap", "RColorBrewer",
  "circlize", "ragg", "pROC", "RANN", "Rtsne", "enrichR", "BiocManager"))
BiocManager::install(c("AUCell", "ComplexHeatmap", "clusterProfiler", "GSVA",
  "limma", "MAST", "scRepertoire"), ask = FALSE, update = FALSE)
```

Check the selected script's `library()` calls for additional packages. CellChat and BAM workflows require separate dependencies. Figure exports use Arial.

## Data and Paths

Extract the required archives beside `Code/`; each ZIP already contains the `Data/` prefix.

```text
project/
|-- Code/
`-- Data/
```

| Data path | Main contents | Download archive |
| --- | --- | --- |
| `Data/synthetic/`, `Data/pbmc30/`, `Data/pbmc4/`, `Data/hash10x/` | Benchmark source tables, HLA typing/coverage, `4mix_PBMC.rds` and `10X_hash_cells.rds`. | [Data_Benchmarks.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/Benchmark_data/Data_Benchmarks.zip) |
| `Data/csf/` | `CSF_scRNA.rds`, `CSF_Myeloid_revision.rds`, `CSF_Tcell.rds`, sample information and plotting sources. | [Data_CSF.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/CSF_data/Data_CSF.zip) |
| `Data/balf/` | `lavage_Tcell.rds`, sample information, CAP/LC signature genes/scores and plotting sources. | [Data_BALF.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/BALF_data/Data_BALF.zip) |
| `Data/urine/` | `Urine_ALL.rds`, `Urine_Tcell.rds`, `Tcell_only_ABMR_TCMR_reclustered.rds`, sample information, signatures and source tables. | [Data_Urine.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/Urine_data/Data_Urine.zip) |
| `Data/gse189536/` | `Immune_annotated.rds` and `Tcell_harmony.rds` for kidney-biopsy comparisons. | [Data_GSE189536.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/Kidney_biopy_data/Data_GSE189536.zip) |
| `Data/references/` | `gene_modules.csv`, `disease_signatures.csv` and GO biological-process gene sets. | [Data_References.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/Reference_data/Data_References.zip) |

Each RDS has a matching `*_metadata.csv.gz`; `sample_metadata.csv` records object-specific sample counts. Preserve file names and capitalization. See [data documentation](../Data/README.md) for details.

## Manual Use

1. Select the script using the relevant Figure-directory README and install its required packages.
2. Edit input, output and `source()` paths inside the script. If it uses `01_setup_and_helpers.R`, update that file's paths as well.
3. Run the required sections manually in R/RStudio, starting with package loading and path settings. Reuse saved objects for downstream figures; script numbers are not an execution order.
4. Check the configured output folders for PDF/PNG figures and tables. Review export paths to avoid overwriting supplied objects or signature tables.

Raw matrices, BAM/VCF files, TCR contigs and signature-search intermediates are not included; upstream reconstruction requires additional inputs. There is no one-command runner.

For BAM benchmarks, edit `Figure1_benchmark/BAM_demultiplex/config.sh` and follow the [benchmark README](Figure1_benchmark/BAM_demultiplex/README.md).
