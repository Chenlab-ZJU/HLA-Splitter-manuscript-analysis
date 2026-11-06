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

## Data Sources and Contents

This bundle contains 11 saved Seurat objects and supporting metadata, HLA genotypes, gene sets and figure source tables for downstream analysis.

The clinical cohorts below were generated in this study. Packaged objects represent the analyzed cells/subsets, not every collected specimen.

| Dataset | Data source | Data type and main contents | Download archive |
| --- | --- | --- | --- |
| Synthetic benchmark | Synthetic pools constructed from 40 individual PBMC scRNA-seq datasets in [GSE226602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226602), an Alzheimer's disease peripheral-immunity resource. | Benchmark result tables: `Cellnum_benchmark.csv`, `Demultiplex_performance.csv`, `SeqDepth_Benchmark.csv`, HLA-resolution/restriction comparisons, HLA abundance, cell-level assignments and donor genotypes. | [Data_Benchmarks.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/Benchmark_data/Data_Benchmarks.zip) |
| PBMC30 | Experimental 30-donor PBMC pool generated in this study. | `HLA_typing.csv`, `Genome_coverage.csv`, `HLA_locus_coverage.csv` and the saved ONT/SBT HLA embedding/assignment source. No PBMC30 expression object is included. | Included in the benchmark archive |
| PBMC4 | Experimental four-donor sorted-cell pool generated in this study, with T cells, NK cells, monocytes and B cells from different donors. | `4mix_PBMC.rds`, donor HLA genotypes and cell/sample metadata. | Included in the benchmark archive |
| 10x hashing | Public [10x Genomics 5' immune-profiling hashing dataset](https://www.10xgenomics.com/datasets/5-hashing-example-with-tabs-2-standard), comprising four donors and PBMC/BMMC samples. | `10X_hash_cells.rds`, donor HLA genotypes and cell/sample metadata for comparing HLA assignments with hashtag labels. | Included in the benchmark archive |
| CSF | Clinical cerebrospinal-fluid samples collected at The Second Affiliated Hospital of Zhejiang University School of Medicine in this study; metastases, hemorrhage and encephalitis groups. | `CSF_scRNA.rds`, `CSF_Myeloid_revision.rds`, `CSF_Tcell.rds`, coded sample information, HLA genotypes, composition/module/key-gene source tables and the manual T-cell annotation map. | [Data_CSF.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/CSF_data/Data_CSF.zip) |
| BALF | Clinical bronchoalveolar-lavage samples collected at The Second Affiliated Hospital of Zhejiang University School of Medicine in this study; downstream T-cell comparisons use CAP, LC and ILD groups. | `lavage_Tcell.rds`, `SampleInfo.csv`, `Geneset.csv`, `Sample_Score.csv`, subtype DEG/color sources, HLA assignments and saved embedding source. | [Data_BALF.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/BALF_data/Data_BALF.zip) |
| Urine | Clinical urine samples from kidney-transplant recipients at The First Affiliated Hospital of Zhejiang University in this study; NR, TCMR, ABMR, mixed rejection and BK groups. | `Urine_ALL.rds`, `Urine_Tcell.rds`, `Tcell_only_ABMR_TCMR_reclustered.rds`, coded sample information, HLA genotypes, DEG/signature tables and demultiplexing/immune-cell-recovery sources. | [Data_Urine.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/Urine_data/Data_Urine.zip) |
| GSE189536 | Public [kidney-allograft biopsy CITE-seq dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189536) from Leckie-Harre et al., [Transplantation (2024)](https://pubmed.ncbi.nlm.nih.gov/37638864/), covering 11 biopsies with non-rejection, TCMR, ABMR or mixed rejection. | Reprocessed RNA-analysis objects `Immune_annotated.rds` and `Tcell_harmony.rds`, plus cell/sample metadata. These are this study's processed immune/T-cell subsets, not the full original CITE-seq deposit. | [Data_GSE189536.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/Kidney_biopy_data/Data_GSE189536.zip) |
| References | Functional modules curated in the analysis scripts, final disease signatures selected in this study, and the [MSigDB C5 GO biological-process collection](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C5). | `gene_modules.csv`, `disease_signatures.csv` and `c5.go.bp.v2023.2.Hs.symbols.gmt` (human gene symbols, release 2023.2). | [Data_References.zip](https://github.com/Chenlab-ZJU/HLA-Splitter-manuscript-analysis/releases/download/Reference_data/Data_References.zip) |

CAP: community-acquired pneumonia; LC: lung cancer; ILD: interstitial lung disease; NR: non-rejection; TCMR: T cell-mediated rejection; ABMR: antibody-mediated rejection; BK: BK polyomavirus nephropathy.

## Data and Paths

<details>
<summary><code>project/</code></summary>
<blockquote>
<details>
<summary><code>Code/</code></summary>
<blockquote>

```text
DATA_README.md
README.md
```

<details>
<summary><code>00_Shared/</code></summary>
<blockquote>

```text
enrichment_analysis.R
export_minimal_data.R
plot_cell_proportion_barplot.R
plot_cell_ratio_comparison.R
plot_celltype_ratio.R
plot_dotplot.R
plot_grouped_umap.R
plot_hla_tsne.R
plot_pseudobulk_heatmap.R
plot_top_pseudobulk.R
plot_violin_comparison.R
plot_volcano.R
```

</blockquote>
</details>
<details>
<summary><code>Figure1_benchmark/</code></summary>
<blockquote>

```text
01_hla_demultiplex.R
02_demultiplex_evaluation.R
03_synthetic_scrna_preprocessing.R
04_synthetic_pool_evaluation.R
05_cell_number_depth_hla_restriction.R
06_tool_accuracy_recovery_runtime.R
07_hla_expression_and_accuracy.R
08_hla_gene_abundance.R
09_nanopore_hla_typing.R
10_pbmc30_validation.R
11_hla_resolution_comparison.R
12_nanopore_hla_accuracy.R
13_nanopore_hla_coverage.R
14_nanopore_sbt_assignment.R
15_pbmc4_validation.R
16_hashing_validation.R
README.md
```

<details>
<summary><code>BAM_demultiplex/</code></summary>
<blockquote>

```text
README.md
config.sh
```

<details>
<summary><code>Benchmarking/</code></summary>
<blockquote>

```text
03_run_benchmark.sh
```

</blockquote>
</details>
<details>
<summary><code>Data_Simulation/</code></summary>
<blockquote>

```text
01_simulate_pools.sh
```

</blockquote>
</details>
<details>
<summary><code>Evaluation/</code></summary>
<blockquote>

```text
04_evaluate_predictions.py
```

</blockquote>
</details>
<details>
<summary><code>Reference_Preparation/</code></summary>
<blockquote>

```text
02_prepare_references.sh
bam_to_vcf.py
extract_hla_genotypes.py
```

</blockquote>
</details>
<details>
<summary><code>environment_setup/</code></summary>
<blockquote>

```text
setup_environment.sh
```

</blockquote>
</details>
</blockquote>
</details>
</blockquote>
</details>
<details>
<summary><code>Figure2_CSF/</code></summary>
<blockquote>

```text
01_setup_and_helpers.R
02_data_preprocessing.R
03_demultiplexing_and_doublets.R
04_global_clustering.R
05_myeloid_analysis.R
06_tcell_analysis.R
07_scrna_celltype_composition_plots.R
08_myeloid_reanalysis.R
09_myeloid_group_deg_function.R
10_tcell_subtype_function_tcr.R
11_scrna_all_celltype_facet.R
12_selected_celltype_proportions.R
13_myeloid_annotation_umap.R
14_myeloid_gene_module_heatmap.R
15_myeloid_module_scores.R
16_scrna_all_major_celltype_heatmap.R
17_myeloid_gene_expression.R
18_myeloid_key_gene_feature_umap.R
19_tcell_annotation_umap.R
20_tcell_celltype_proportion_heatmap.R
21_tcell_marker_heatmap.R
22_tcell_subtype_proportion_compare.R
23_tcell_gene_expression.R
24_tcell_functional_module_score_compare.R
README.md
```

</blockquote>
</details>
<details>
<summary><code>Figure3_BALF/</code></summary>
<blockquote>

```text
01_setup_and_helpers.R
02_data_preprocessing_demux.R
03_demultiplexing_and_doublets.R
04_global_clustering.R
05_tcr_repertoire.R
06_aucell_scoring.R
07_pseudobulk_functional_analysis.R
08_disease_signature_selection.R
09_tcell_group_umap_composition.R
10_tcell_subtype_proportions.R
11_tcell_subtype_deg_heatmap.R
12_tcell_functional_module_scores.R
13_signature_scores_roc.R
14_tcell_module_umap.R
15_hla_assignment_overview.R
16_cap_lc_tcell_deg.R
17_tcr_expansion_summary.R
README.md
```

</blockquote>
</details>
<details>
<summary><code>Figure4_urine_immune/</code></summary>
<blockquote>

```text
01_setup_and_helpers.R
02_data_integration_and_qc.R
03_demultiplexing.R
04_global_clustering.R
05_aucell_scoring.R
06_pseudobulk_analysis.R
07_celltype_umap_and_rejection_proportions.R
08_major_celltype_group_umap.R
09_group_celltype_composition.R
10_demultiplex_immune_cell_ratio_plots.R
11_all_cells_rejection_module_heatmap_score.R
12_demultiplex_hla_tsne_predict_sample_number.R
13_sample_celltype_proportions.R
14_major_celltype_hla_cytotoxic_exhaustion_violin.R
15_signature_score_analysis_gene_contribution.R
16_signature_candidate_evaluation.R
17_signature_cross_cohort_robust_selection.R
18_cross_cohort_signature_distribution.R
19_cross_cohort_gene_correlation.R
20_cross_cohort_group_comparisons.R
README.md
```

</blockquote>
</details>
<details>
<summary><code>Figure5_urine_Tcell_CXCR3/</code></summary>
<blockquote>

```text
01_tcell_subclustering.R
02_tcr_repertoire_analysis.R
03_cell_communication.R
04_tcell_recluster_subtype_frequency.R
05_abmr_vs_tcmr_tcell_de_function.R
06_balanced_cxcr3_ligand_communication.R
07_tcell_group_umap.R
08_tcell_group_deg_module_heatmap_score.R
09_cxcr3_expression_positive_fraction.R
10_subtype_cxcr3_positive_fraction.R
11_expanded_tcell_fraction_by_subtype.R
12_cxcr3_expanded_tcell_fraction.R
13_tcr_expansion_umap_distribution.R
14_expanded_tcr_clone_diversity_by_group.R
15_gse189536_immune_annotation.R
16_gse189536_rejection_axis_validation.R
17_gse189536_cxcr3_marker_dotplot.R
18_gse189536_tcell_modules.R
README.md
```

</blockquote>
</details>
</blockquote>
</details>
<details>
<summary><code>Data/</code></summary>
<blockquote>
<details>
<summary><code>balf/</code></summary>
<blockquote>

```text
Geneset.csv
HLA-genotypes.csv
HLA_TSNE_Adjusted_Prediction_revision_source.csv
HLA_demultiplex.csv
SampleInfo.csv
Sample_Score.csv
lavage_Tcell.rds
lavage_Tcell_UMAP_Cluster_Group_color_map.csv
lavage_Tcell_metadata.csv.gz
lavage_Tcell_subtype_DEG_wilcox_revision.csv
sample_metadata.csv
```

</blockquote>
</details>
<details>
<summary><code>csf/</code></summary>
<blockquote>

```text
CSF_Myeloid_revision.rds
CSF_Myeloid_revision_metadata.csv.gz
CSF_Tcell.rds
CSF_Tcell_metadata.csv.gz
CSF_scRNA.rds
CSF_scRNA_metadata.csv.gz
HLA-genotypes.csv
Myeloid_revision_functional_module_sample_scores.csv
Myeloid_revision_group_DEG_heatmap_source.csv
Myeloid_revision_key_gene_sample_expression.csv
SampleInfo.csv
Tcell_revision_manual_cluster_reannotation.csv
sample_metadata.csv
scRNA_celltype_sample_percent_revision.csv
scRNA_major_celltype_sample_percent_revision.csv
```

</blockquote>
</details>
<details>
<summary><code>gse189536/</code></summary>
<blockquote>

```text
Immune_annotated.rds
Immune_annotated_metadata.csv.gz
Tcell_harmony.rds
Tcell_harmony_metadata.csv.gz
sample_metadata.csv
```

</blockquote>
</details>
<details>
<summary><code>hash10x/</code></summary>
<blockquote>

```text
10X_hash_cells.rds
10X_hash_cells_metadata.csv.gz
HLA-genotypes.csv
sample_metadata.csv
```

</blockquote>
</details>
<details>
<summary><code>pbmc30/</code></summary>
<blockquote>

```text
Genome_coverage.csv
HLA_TSNE_ONT_SBT_revision_source.csv
HLA_locus_coverage.csv
HLA_typing.csv
```

</blockquote>
</details>
<details>
<summary><code>pbmc4/</code></summary>
<blockquote>

```text
4mix_PBMC.rds
4mix_PBMC_metadata.csv.gz
HLA-genotypes.csv
sample_metadata.csv
```

</blockquote>
</details>
<details>
<summary><code>references/</code></summary>
<blockquote>

```text
c5.go.bp.v2023.2.Hs.symbols.gmt
disease_signatures.csv
gene_modules.csv
```

</blockquote>
</details>
<details>
<summary><code>synthetic/</code></summary>
<blockquote>

```text
Cellnum_benchmark.csv
Demultiplex_performance.csv
Figure_S1e_HLA_gene_abundance_source.csv
HLA-genotypes.csv
HLA_Resolution_benchmark.csv
HLA_Restricted_benchmark.csv
SeqDepth_Benchmark.csv
demultiplex_result.csv
```

</blockquote>
</details>
<details>
<summary><code>urine/</code></summary>
<blockquote>

```text
Comparision_immune_cell_ratio.csv
Cross_cohort_robust_augmented_rejection_gene_sets_used_revision.csv
Cross_cohort_robust_sample_scores.csv
Cross_cohort_robust_selected_genes.csv
HLA-genotypes.csv
HLA_TSNE_Adjusted_Prediction_source.csv
Predict_sample_number.csv
Sample_Info.csv
Signature_gene_sets_used_revision.csv
Tcell_Group_DEG.csv
Tcell_only_ABMR_TCMR_reclustered.rds
Tcell_only_ABMR_TCMR_reclustered_metadata.csv.gz
Urine_ALL.rds
Urine_ALL_metadata.csv.gz
Urine_Tcell.rds
Urine_Tcell_metadata.csv.gz
sample_metadata.csv
scRNA_ABMR_DEG.csv
scRNA_BK_DEG.csv
scRNA_TCMR_DEG.csv
```

</blockquote>
</details>
</blockquote>
</details>
</blockquote>
</details>

## Manual Use

1. Select the script using the relevant Figure-directory README and install its required packages.
2. Edit input, output and `source()` paths inside the script. If it uses `01_setup_and_helpers.R`, update that file's paths as well.
3. Run the required sections manually in R/RStudio, starting with package loading and path settings. Reuse saved objects for downstream figures; script numbers are not an execution order.
4. Check the configured output folders for PDF/PNG figures and tables. Review export paths to avoid overwriting supplied objects or signature tables.

Raw matrices, BAM/VCF files, TCR contigs and signature-search intermediates are not included; upstream reconstruction requires additional inputs. There is no one-command runner.

For BAM benchmarks, edit `Figure1_benchmark/BAM_demultiplex/config.sh` and follow the [benchmark README](Figure1_benchmark/BAM_demultiplex/README.md).
