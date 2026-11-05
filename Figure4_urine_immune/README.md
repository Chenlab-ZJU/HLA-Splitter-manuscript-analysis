# Figure4 urine immune

Also supports Supplementary Figures S7-S8.

Open the required R script, edit its input/output and helper-script paths, and run the required sections manually. Upstream scripts require original inputs and are not needed to redraw figures from saved objects. Numbers are an index, not an automatic execution order.

| Script | Role |
| --- | --- |
| `01_setup_and_helpers.R` | functions: setup and helpers |
| `02_data_integration_and_qc.R` | upstream: data integration and QC |
| `03_demultiplexing.R` | upstream: demultiplexing |
| `04_global_clustering.R` | upstream: global clustering |
| `05_aucell_scoring.R` | upstream: AUCell scoring |
| `06_pseudobulk_analysis.R` | upstream: Pseudobulk analysis |
| `07_celltype_umap_and_rejection_proportions.R` | analysis_or_plot: celltype UMAP and rejection proportions |
| `08_major_celltype_group_umap.R` | analysis_or_plot: major celltype group UMAP |
| `09_group_celltype_composition.R` | analysis_or_plot: group celltype composition |
| `10_demultiplex_immune_cell_ratio_plots.R` | analysis_or_plot: demultiplex immune cell ratio plots |
| `11_all_cells_rejection_module_heatmap_score.R` | analysis_or_plot: all cells rejection module heatmap score |
| `12_demultiplex_hla_tsne_predict_sample_number.R` | analysis_or_plot: demultiplex HLA tsne predict sample number |
| `13_sample_celltype_proportions.R` | analysis_or_plot: sample celltype proportions |
| `14_major_celltype_hla_cytotoxic_exhaustion_violin.R` | analysis_or_plot: major celltype HLA cytotoxic exhaustion violin |
| `15_signature_score_analysis_gene_contribution.R` | analysis_or_plot: signature score analysis gene contribution |
| `16_signature_candidate_evaluation.R` | analysis_or_plot: signature score advanced optimization GSE189536 validation |
| `17_signature_cross_cohort_robust_selection.R` | analysis_or_plot: signature cross cohort robust selection |
| `18_cross_cohort_signature_distribution.R` | analysis_or_plot: cross cohort robust signature distribution plots |
| `19_cross_cohort_gene_correlation.R` | analysis_or_plot: cross cohort signature gene expression correlation |
| `20_cross_cohort_group_comparisons.R` | analysis_or_plot: cross cohort robust boxplots pairwise target vs each group |
