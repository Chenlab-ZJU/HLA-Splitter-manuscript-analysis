# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 03_global_clustering.R
# Data cohort: balf. Inputs use the compact English data layout.

# ==============================================================================
# File: 02_global_clustering.R
# Description: global clustering, cell-type annotation and T-cell subclustering.
# ==============================================================================
rm(list = ls())
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure3_BALF/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA_Intergrate = readRDS("Result/lavage_scRNA_Intergrate.rds")
#Tcell = readRDS("./Result/lavage_Tcell.rds")

# Select cells and recluster.
scRNA <- subset(scRNA_Intergrate, subset = !HLA_Splitter %in% c("Doublets", "Undefined"))
SampleInfo <- read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "balf", "SampleInfo.csv"))
Sample_levels <- c(mixedsort(SampleInfo$Patient.ID))
scRNA$Group = factor(scRNA$Group, setdiff(scRNA$Group,c("Undefined","Doublet")))
# Exclude samples with fewer than 50 cells.
valid_samples <- Sample_levels[table(scRNA$HLA_Splitter) > 50]
scRNA <- subset(scRNA, subset = HLA_Splitter %in% valid_samples)

# Cluster and annotate the global cell population.
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>%
  RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.8)

# Identify and visualize global markers.
export_top_markers(scRNA, "seurat_clusters", "scRNA")

# Manual checkpoint: fill the Celltype column in Result/scRNA_makers_summary.csv before continuing.
df <- read.csv("./Result/scRNA_makers_summary.csv")
Celltype <- df$Celltype
names(Celltype) <- as.character(df$seurat_clusters)
scRNA <- SetIdent(scRNA, value = "seurat_clusters")
scRNA <- RenameIdents(scRNA, Celltype)
scRNA$Celltype <- Idents(scRNA)
unique(scRNA$Celltype) %>% print()

# This section assumes that Celltype annotations have already been assigned.
Tcell = scRNA[,scRNA$Celltype != "Contament"]
Tcell$Celltype <- factor(Tcell$Celltype, levels = c("Naive T", "Quiescent T", "CD4+ Tact", "Th17-like", "Tph", "IFN+ T", "Tregs", 
                                                    "CD8+ Tem", "CD8+ Trm", "CD8+ Teff", "Cycling T", "MAIT cells"))

# Recluster the selected cells.
Tcell = Tcell[,Tcell$Group != "Others"]
Tcell$Group = factor(Tcell$Group, levels = c("CAP","LC","ILD"))
                       
Tcell <- RunUMAP(Tcell,dims = 1:30,  min.dist = 0.1)
#FeaturePlot(Tcell, features = c("CD4","CD8A","KLRD1","ISG15"))

# Draw the three-panel global UMAP.
p1 <- DimPlot(Tcell, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(Tcell, group.by = "Celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p3 <- DimPlot(Tcell, group.by = "Group", label = TRUE, alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
ggsave("./Figures/Tcell_UMAP_celltype.pdf", p1 + p2 + p3, height = 5, width = 18)

# Plot global marker genes.
global_mk <- c("CCR7", "LEF1", "BRD1", "GRK5", "CD40LG", "HLA-DR", 
               "CCR6", "KLRB1", "CXCL13", "TNFRSF18", "IFIT3", "ISG15", "FOXP3", "IL2RA",
               "GZMK", "EOMES", "ZNF683","ITGAE",  "GZMH", "PRF1", 
               "PDCD1", "MKI67","TYMS", "SLC4A10", "TRAV1-2")
p_dot <- DotPlot(Tcell, group.by = "Celltype", features = global_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x=NULL, y=NULL) +
  scale_color_gradientn(colours = c("navy","white","firebrick3"))
ggsave("./Figures/Tcell_Dotplot_clusters.pdf", plot = p_dot, height = 4, width = 9)

# Draw group-faceted UMAPs.
combined_plot <- ggarrange(plotlist = plot_grouped_umap(Tcell, "Group"), ncol = 3,common.legend = TRUE, legend = "right")

ggsave("./Figures/Tcell_all_celltype_facet.pdf", combined_plot, height = 5, width = 12)

# Draw the cell-type heatmap.
p_list <- plot_celltype_ratio(seurat_obj = Tcell,
                              sample_col = "HLA_Splitter", group_col = "Group", cluster_col = "Celltype",
                              comparison_groups = comparison, group_colors = GROUP_COLORS, cluster_colors = GLOBAL_COLORS)

ggsave("./Figures/Tcell_celltype_stacked_bar.pdf", p_list$stacked_bar, width = 5, height = 4)
ggsave("./Figures/Tcell_celltype_heatmap.pdf", p_list$heatmap, width = 5, height = 5)


# Plot cell-subtype proportions.
comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))
p_cell_ratio <- plot_cell_ratio_comparison(cell_subset = Tcell,
                                           sample_col = "HLA_Splitter", cluster_col =  "Celltype", group_col = "Group",
                                           facet = T, comparison_groups = comparison,
                                           write_to_csv = "./Result/Tcell_celltype_ratio", plot_colors = GLOBAL_COLORS)
ggsave("./Figures/Tcell_celltype_compare.pdf", p_cell_ratio, width = 10, height = 10)

# Save analysis outputs.
# Export the analysis source data for figure reproduction.
saveRDS(Tcell, file.path("D:/项目/2026.Nature Communications/Revison/Data", "balf", "lavage_Tcell.rds"))

# ==============================================================================
# Description: compare cell composition, differential expression and enrichment across groups.
# ==============================================================================
Tcell = readRDS(file.path("D:/项目/2026.Nature Communications/Revison/Data", "balf", "lavage_Tcell.rds"))
Tcell = Tcell[,Tcell$Group!= "Others"]

comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))

# ==============================================================================
# Export pseudobulk lollipop plots, group-specific CSVs and bar/dot plots.
# ==============================================================================
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_top_pseudobulk.R"))

out_dir = "./Result/Marker_Analysis"

plot_top_pseudobulk(Tcell, 
                    group_by = "Group",
                    sample_col = "HLA_Splitter",
                    top_n_lollipop = 10, 
                    top_n_bardot = 10, 
                    out_dir = out_dir,
                    my_colors = GROUP_COLORS)

unique_groups = unique(Tcell$Group)

Group_DEG = data.frame()
for (grp in unique_groups) {
  top_gene_tmp = read.csv(paste0(out_dir,"/Marker_Stats_", grp, ".csv"))
  Group_DEG = rbind(Group_DEG, top_gene_tmp)
}
write.csv(Group_DEG, paste0(out_dir, "/Tcell_Group_DEG.csv"), row.names = FALSE)

# ==============================================================================
# Combine top group-associated genes in a pseudobulk heatmap.
# ==============================================================================
message("Generating global pseudobulk heatmap...")

# Select genes with the highest log2FC in each group.
top_heatmap_genes <- Group_DEG %>%
  group_by(Group_Col) %>%
  dplyr::slice_max(order_by = PB_Log2FC, n = 5, with_ties = FALSE) %>%
  pull(gene) %>%
  unique()

Tcell_avg <- AverageExpression(Tcell, group.by = "Group", assays = "RNA", return.seurat = TRUE)
Tcell_avg <- ScaleData(Tcell_avg, features = rownames(Tcell_avg))
gene_matrix_t <- Tcell_avg@assays$RNA$scale.data %>% as.matrix()

# Keep genes available in the matrix used for the heatmap.
genes_to_plot <- intersect(top_heatmap_genes, rownames(gene_matrix_t))

plot_heatmap_t_bulk <- pheatmap(gene_matrix_t[genes_to_plot, ], 
                                cluster_rows = TRUE, cluster_cols = FALSE, 
                                angle_col = 45, cellwidth = 40, cellheight = 12) # Adjust cellheight to limit the heatmap length.

ggsave("./Figures/Tcell_Heatmap_bulk_AllGroups.pdf", plot_heatmap_t_bulk, height = 12, width = 8)

# ==============================================================================
# Draw group-specific violin plots and enrichment results.
# ==============================================================================
unique_groups <- as.character(unique(Group_DEG$Group_Col))

for (grp in unique_groups) {
  message("\n>> Analyzing subset: ", grp)
  
  # Draw violin plots.
  plot_top_violins(
    seu_obj = Tcell,
    deg_data = Group_DEG,
    target_group = grp,
    out_dir = out_dir,
    comparisons_list = comparison,
    plot_colors = GROUP_COLORS,
    group_col = "Group_Col",
    lfc_col = "PB_Log2FC",
    lfc_threshold = 0.1,
    top_n = 8
  )
  
  # Run GO enrichment and plot results.
  run_go_enrichment(
    deg_data = Group_DEG,
    target_group = grp,
    out_dir = out_dir,
    group_col = "Group_Col",
    lfc_col = "PB_Log2FC",
    lfc_threshold = 0.1,
    top_n_pathways = 8
  )
}

# ==============================================================================
# Run the combined pseudobulk differential-expression analysis.
# ==============================================================================
comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))

source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure3_BALF/07_pseudobulk_functional_analysis.R"))

#Tcell$Major_type = "Tcell"

res <- run_sc_pseudobulk_functional_analysis(
  seu = Tcell,
  group_col = "Group",
  sample_col = "HLA_Splitter",
  celltype_col = "Major_type",
  default_celltype = "Tcell",
  comparisons = comparison,
  organism = "human",
  gene_id_type = "SYMBOL",
  assay = "RNA",
  layer = "counts",
  min_cells_per_sample = 20,
  min_samples_per_group = 2,
  min_count = 20,
  min_samples_expr = 3,
  padj_cutoff = 0.05,
  lfc_cutoff = 0.25,
  output_dir = "./Result/Tcell_pseudobulk_DEG_GSEA"
)

# Draw a heatmap using the specified p-value criterion.
heatmap_plot <- plot_lrt_global_heatmap(
  lrt_table = res$overall_LRT,
  pb_counts = res$pseudobulk$counts,
  pb_meta = res$pseudobulk$metadata,
  target_celltype = "Tcell",
  filter_by = "pvalue",    # Use unadjusted p values explicitly.
  p_cutoff = 0.01,         # Use the specified 0.01 p-value cutoff.
  k_clusters = 4
)

# Cluster expression trends using the specified p-value criterion.
trend_res <- plot_lrt_cluster_trends(
  lrt_table = res$overall_LRT,
  pb_counts = res$pseudobulk$counts,
  pb_meta = res$pseudobulk$metadata,
  target_celltype = "Tcell",
  filter_by = "pvalue",    # Use unadjusted p values explicitly.
  p_cutoff = 0.01,         # Use the specified 0.01 p-value cutoff.
  k_clusters = 4
)

# Display the resulting plot.
print(trend_res$plot)
