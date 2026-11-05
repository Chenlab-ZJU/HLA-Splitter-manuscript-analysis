# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream cohort processing and reusable helper definitions.
# Original source: 00_setup_and_helpers.R
# Data cohort: csf. Inputs use the compact English data layout.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.

# ==============================================================================
# File: 00_setup_and_helpers.R
# Description: load dependencies, configure paths and source helper functions.
# ==============================================================================

WORK_DIR <- file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf")
setwd(WORK_DIR)
set.seed(2026) # Set the random seed for reproducible stochastic analysis steps.

# Load required packages.
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(harmony)
  library(reshape2)
  library(pheatmap)
  library(sctransform)
  library(RColorBrewer)
  library(grid)
  library(stringr)
  library(gtools)
  library(ggsci)
  library(ggpubr)
  library(patchwork)
  library(scrubletR)
  library(paletteer)
  library(ggExtra)
  library(enrichR)
  library(GSVA)
  library(msigdbr)
  library(clusterProfiler)
  library(limma)
  library(ggrepel)
  library(presto)
})

# Configure named color palettes.
GLOBAL_COLORS <- c(
  "#E42320","#699ECA","#4DAF4A","#FF8C00","#F898CB","#D65190","#731A73","#FFCB5B",
  "#E89B1E","#0076B9","#3D505A","#0098B2","#8A7355","#CC5B45","#008B8B","#4E5689",
  "#F9B3AD","#EDC66A","#9FDAF7","#C9A1CA","#D9DEE7","#C74546","#EB7E60","#8FB4DC"
)
GROUP_COLORS <- c("Metastases"= "#E42320", "Hemorrhage" = "#699ECA", "Encephalitis" = "#4DAF4A", "Doublets" = "black", "Undefined" = "grey")

comparison <- list(c("Metastases","Hemorrhage"), c("Metastases","Encephalitis"), c("Hemorrhage","Encephalitis"))

# Source custom helper scripts.
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure1_benchmark/01_hla_demultiplex.R"))
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_grouped_umap.R"))
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_cell_proportion_barplot.R"))
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_hla_tsne.R"))
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_celltype_ratio.R"))
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_cell_ratio_comparison.R"))
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_volcano.R"))
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/enrichment_analysis.R"))
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_pseudobulk_heatmap.R"))

# Define helper functions.
# Export the top 10 subtype differential-expression markers.
export_top_markers <- function(seurat_obj, ident_col, out_prefix) {
  Idents(seurat_obj) <- ident_col
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Export the analysis source data for figure reproduction.
  write.csv(markers, paste0("./Result/", out_prefix, "_makers.csv"), row.names = FALSE)
  
  df_summary <- seurat_obj@meta.data %>%
    select(all_of(ident_col)) %>% distinct() %>% arrange(across(all_of(ident_col)))
  
  df_top_genes <- markers %>%
    filter(p_val_adj < 0.01) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>%
    summarise(top_genes = paste0(gene, collapse = ", ")) %>%
    rename(!!ident_col := cluster)
  
  final_df <- left_join(df_summary, df_top_genes, by = ident_col)
  final_df$Celltype = c("")
  write.csv(final_df, paste0("./Result/", out_prefix, "_makers_summary.csv"), row.names = FALSE)
  
  message("Open ./Result/", out_prefix, "_makers_summary.csv, add the Celltype column, and save before continuing.")
}
