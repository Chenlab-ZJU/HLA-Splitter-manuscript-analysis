# Figure2 CSF

Also supports Supplementary Figures S4-S5.

Open the required R script, edit its input/output and helper-script paths, and run the required sections manually. Upstream scripts require original inputs and are not needed to redraw figures from saved objects. Numbers are an index, not an automatic execution order.

| Script | Role |
| --- | --- |
| `01_setup_and_helpers.R` | functions: setup and helpers |
| `02_data_preprocessing.R` | upstream: data preprocessing |
| `03_demultiplexing_and_doublets.R` | upstream: demultiplexing and doublets |
| `04_global_clustering.R` | upstream: global clustering |
| `05_myeloid_analysis.R` | upstream: myeloid analysis |
| `06_tcell_analysis.R` | upstream: tcell analysis |
| `07_scrna_celltype_composition_plots.R` | analysis: scRNA celltype composition plots |
| `08_myeloid_reanalysis.R` | analysis: myeloid reanalysis |
| `09_myeloid_group_deg_function.R` | analysis: myeloid group DEG function |
| `10_tcell_subtype_function_tcr.R` | analysis: tcell subtype function tcr |
| `11_scrna_all_celltype_facet.R` | analysis_or_plot: scRNA all celltype facet |
| `12_selected_celltype_proportions.R` | analysis_or_plot: selected celltype proportions |
| `13_myeloid_annotation_umap.R` | analysis_or_plot: myeloid annotation UMAP |
| `14_myeloid_gene_module_heatmap.R` | analysis_or_plot: myeloid gene module heatmap |
| `15_myeloid_module_scores.R` | analysis_or_plot: myeloid module scores |
| `16_scrna_all_major_celltype_heatmap.R` | analysis_or_plot: scRNA all major celltype heatmap |
| `17_myeloid_gene_expression.R` | analysis_or_plot: myeloid gene expression |
| `18_myeloid_key_gene_feature_umap.R` | analysis_or_plot: myeloid key gene feature umap |
| `19_tcell_annotation_umap.R` | analysis_or_plot: Tcell annotation UMAP |
| `20_tcell_celltype_proportion_heatmap.R` | analysis_or_plot: tcell celltype proportion heatmap |
| `21_tcell_marker_heatmap.R` | analysis_or_plot: Tcell marker heatmap |
| `22_tcell_subtype_proportion_compare.R` | analysis_or_plot: tcell subtype proportion compare |
| `23_tcell_gene_expression.R` | analysis_or_plot: Tcell gene expression |
| `24_tcell_functional_module_score_compare.R` | analysis_or_plot: tcell functional module score compare |
