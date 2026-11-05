# Figure1 benchmark

Upstream BAM-based tool execution is in [BAM_demultiplex](BAM_demultiplex/README.md). Run those stages manually only when new benchmarking is required; they are not needed to redraw figures from saved source tables.

Also supports Supplementary Figures S1-S3.

Open the required R script, edit its input/output and helper-script paths, and run the required sections manually. Upstream scripts require original inputs and are not needed to redraw figures from saved objects. Numbers are an index, not an automatic execution order.

| Script | Role |
| --- | --- |
| `01_hla_demultiplex.R` | functions: HLA demultiplex v2 |
| `02_demultiplex_evaluation.R` | functions: Demultiplex evaluation |
| `03_synthetic_scrna_preprocessing.R` | upstream: scRNA analysis |
| `04_synthetic_pool_evaluation.R` | upstream: Syn data |
| `05_cell_number_depth_hla_restriction.R` | analysis_or_plot: cell number depth HLA restriction |
| `06_tool_accuracy_recovery_runtime.R` | analysis_or_plot: tool accuracy recovery runtime |
| `07_hla_expression_and_accuracy.R` | analysis_or_plot: HLA expression and accuracy |
| `08_hla_gene_abundance.R` | analysis_or_plot: HLA gene abundance |
| `09_nanopore_hla_typing.R` | upstream: Nanopore |
| `10_pbmc30_validation.R` | upstream: 30PBMC |
| `11_hla_resolution_comparison.R` | upstream: Compare HLA resolution |
| `12_nanopore_hla_accuracy.R` | analysis_or_plot: nanopore HLA accuracy |
| `13_nanopore_hla_coverage.R` | analysis_or_plot: nanopore hla locus coverage violin result |
| `14_nanopore_sbt_assignment.R` | analysis_or_plot: nanopore SBT assignment |
| `15_pbmc4_validation.R` | analysis_or_plot: 4mix PBMC |
| `16_hashing_validation.R` | analysis_or_plot: 10X-hash |
