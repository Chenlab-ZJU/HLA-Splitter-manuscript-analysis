# Figure5 urine Tcell CXCR3

Also supports Supplementary Figures S9-S10.

Open the required R script, edit its input/output and helper-script paths, and run the required sections manually. Upstream scripts require original inputs and are not needed to redraw figures from saved objects. Numbers are an index, not an automatic execution order.

| Script | Role |
| --- | --- |
| `01_tcell_subclustering.R` | upstream: Tcell subclustering |
| `02_tcr_repertoire_analysis.R` | upstream: TCR repertoire analysis |
| `03_cell_communication.R` | upstream: Cell communication |
| `04_tcell_recluster_subtype_frequency.R` | upstream: Tcell recluster subtype frequency |
| `05_abmr_vs_tcmr_tcell_de_function.R` | upstream: ABMR vs TCMR Tcell DE function |
| `06_balanced_cxcr3_ligand_communication.R` | upstream: balanced CXCR3 CXCL9 CXCL10 LR communication |
| `07_tcell_group_umap.R` | analysis_or_plot: Tcell group UMAP |
| `08_tcell_group_deg_module_heatmap_score.R` | analysis_or_plot: tcell group DEG module heatmap score |
| `09_cxcr3_expression_positive_fraction.R` | analysis_or_plot: CXCR3 expression positive fraction |
| `10_subtype_cxcr3_positive_fraction.R` | analysis_or_plot: tcell CXCR3 positive ratio by subtype barplot |
| `11_expanded_tcell_fraction_by_subtype.R` | analysis_or_plot: expanded Tcell fraction by subtype |
| `12_cxcr3_expanded_tcell_fraction.R` | analysis_or_plot: CXCR3 expanded Tcell fraction |
| `13_tcr_expansion_umap_distribution.R` | analysis_or_plot: TCR clone ratio UMAP histogram group fraction |
| `14_expanded_tcr_clone_diversity_by_group.R` | analysis_or_plot: expanded TCR clone diversity by group |
| `15_gse189536_immune_annotation.R` | analysis_or_plot: annotate GSE189536 immune like urine |
| `16_gse189536_rejection_axis_validation.R` | analysis_or_plot: GSE189536 rejection axis validation |
| `17_gse189536_cxcr3_marker_dotplot.R` | analysis_or_plot: plot GSE189536 urine style Tcell CXCR3 marker dotplot |
| `18_gse189536_tcell_modules.R` | analysis_or_plot: plot GSE189536 urine style Tcell DEG module heatmap score |
