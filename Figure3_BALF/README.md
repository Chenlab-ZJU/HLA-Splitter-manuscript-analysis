# Figure3 BALF

Also supports Supplementary Figures S6.

Open the required R script, edit its input/output and helper-script paths, and run the required sections manually. Upstream scripts require original inputs and are not needed to redraw figures from saved objects. Numbers are an index, not an automatic execution order.

| Script | Role |
| --- | --- |
| `01_setup_and_helpers.R` | functions: setup and helpers |
| `02_data_preprocessing_demux.R` | upstream: data preprocessing demux |
| `03_demultiplexing_and_doublets.R` | upstream: demultiplexing and doublets |
| `04_global_clustering.R` | upstream: global clustering |
| `05_tcr_repertoire.R` | upstream: tcr repertoire |
| `06_aucell_scoring.R` | upstream: AUCell scoring |
| `07_pseudobulk_functional_analysis.R` | upstream: run sc pseudobulk functional analysis |
| `08_disease_signature_selection.R` | upstream: lavage find signature |
| `09_tcell_group_umap_composition.R` | analysis_or_plot: Tcell group UMAP composition |
| `10_tcell_subtype_proportions.R` | analysis_or_plot: lavage tcell subtype proportion boxplot heatmap |
| `11_tcell_subtype_deg_heatmap.R` | analysis_or_plot: lavage tcell subtype deg module heatmap |
| `12_tcell_functional_module_scores.R` | analysis_or_plot: Tcell functional module scores |
| `13_signature_scores_roc.R` | analysis_or_plot: lavage tcell signature scatter roc contribution |
| `14_tcell_module_umap.R` | analysis_or_plot: Tcell module UMAP |
| `15_hla_assignment_overview.R` | analysis_or_plot: lavage demultiplex hla tsne predict sample number |
| `16_cap_lc_tcell_deg.R` | analysis_or_plot: lavage tcell CAP LC DEG volcano |
| `17_tcr_expansion_summary.R` | analysis_or_plot: lavage tcell tcr clone ratio umap expanded summary |
