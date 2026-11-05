# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream cohort processing and reusable helper definitions.
# Original source: 03_global_clustering.R
# Data cohort: csf. Inputs use the compact English data layout.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.

# ==============================================================================
# File: 03_global_clustering.R
# Description: global dimensional reduction, clustering, annotation and visualization.
# ==============================================================================
rm(list = ls())
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure2_CSF/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA_Intergrate <- readRDS("Result/CSF_scRNA_Intergrate.rds")
#scRNA = readRDS("./Result/CSF_scRNA.rds")

# Select cells and recluster.
scRNA <- subset(scRNA_Intergrate, subset = !HLA_Splitter %in% c("Doublets", "Undefined"))
SampleInfo <- read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "SampleInfo.csv"))
Sample_levels <- c(mixedsort(SampleInfo$Patient.ID))
scRNA$Group = factor(scRNA$Group, setdiff(scRNA$Group,c("Undefined","Doublet")))
# 
# valid_samples <- Sample_levels[table(scRNA$HLA_Splitter) > 100]
# scRNA <- subset(scRNA, subset = HLA_Splitter %in% valid_samples)

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
scRNA = scRNA[,scRNA$Celltype != "Neuronal-like"]
scRNA$Celltype <- factor(scRNA$Celltype, levels = c("Naive T","CD4+ Tact", "IFN+ T", "Tregs", "CD8+ CTL","CD8+ Trm", "CD8+ Tmem",
                                                    "gd T", "B cells", "NK cells","Microglia-like", "Macrophages",
                                                    "Monocytes", "cDC1", "cDC2", "pDC","Proliferating"))

# Draw the three-panel global UMAP.
p1 <- DimPlot(scRNA, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(scRNA, group.by = "Celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p3 <- DimPlot(scRNA, group.by = "Group", label = TRUE, alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
ggsave("./Figures/CSF_UMAP_celltype.pdf", p1 + p2 + p3, height = 5, width = 18)

p_F = FeaturePlot(scRNA, features = c("C1QC", "FCGBP", "TREM2", "SPP1"),ncol=2,cols = c("grey","firebrick"))
ggsave("./Figures/CSF_FeaturePlot.pdf", p_F, height = 8, width = 8)

# Draw group-faceted UMAPs.
combined_plot <- ggarrange(plotlist = plot_grouped_umap(scRNA, "Celltype", "Group", GLOBAL_COLORS), common.legend = TRUE, legend = "right")
ggsave("./Figures/scRNA_all_celltype_facet.pdf", combined_plot, height = 6, width = 8)

# Draw the cell-type heatmap.
p_list <- plot_celltype_ratio(seurat_obj = scRNA,
                              sample_col = "HLA_Splitter", group_col = "Group", cluster_col = "Celltype",
                              comparison_groups = comparison, group_colors = GROUP_COLORS, cluster_colors = GLOBAL_COLORS)

ggsave("./Figures/scRNA_celltype_stacked_bar.pdf", p_list$stacked_bar, width = 5, height = 4)
ggsave("./Figures/scRNA_celltype_heatmap.pdf", p_list$heatmap, width = 5, height = 5)

# Plot global marker genes.
global_mk <- c("TCF7","IL7R", "KLRB1", "CCR6", "IFIT3","ISG15",
               "CTLA4","IL2RA","CD8A", "PRF1", "ZNF683", "FGFBP2","CD27","CISH",
               "XCL1","TRDC","MS4A1","BANK1", "KLRD1", "GNLY", 
               "P2RY12","TREM2", "FCGR3A", "CD163","S100A8","CD14",
               "XCR1", "CLEC9A","CD1C","FCER1A","LILRA4","CLEC4C","MKI67","TYMS")
p_dot <- DotPlot(scRNA, group.by = "Celltype", features = global_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x=NULL, y=NULL) +
  scale_color_gradientn(colours = c("navy","white","firebrick3"))
ggsave("./Figures/CSF_Dotplot_clusters.pdf", plot = p_dot, height = 4, width = 9)

# Plot cell-subtype proportions.
p_cell_ratio <- plot_cell_ratio_comparison(cell_subset = scRNA,
                                           sample_col = "HLA_Splitter", 
                                           cluster_col =  "Celltype",
                                           group_col = "Group",
                                           facet = T, comparison_groups = comparison,
                                           stat_method = "t.test",
                                           write_to_csv = "./Result/scRNA_celltype_ratio", 
                                           plot_colors = GLOBAL_COLORS)
ggsave("./Figures/scRNA_celltype_compare.pdf", p_cell_ratio, width = 10, height = 10)

# Identify differential genes for pairwise disease comparisons.
Idents(scRNA) <- "Group"
# Discover marker genes with the original assay and comparison settings.
CELLDEG <- FindMarkers(scRNA, ident.1 ="Metastases", ident.2 = "Encephalitis", logfc.threshold = 0.1, min.pct = 0.1)
CELLDEG <- CELLDEG[!is.na(CELLDEG$p_val),]
CELLDEG$gene <- rownames(CELLDEG)
# Export the analysis source data for figure reproduction.
write.csv(CELLDEG, "./Result/scRNA_Group_marker_Metastases-Encephalitis.csv")

# Plot differential-expression volcanoes and GO enrichment.
p_volcano <- plot_volcano(diffexp = CELLDEG,
                          top_n_labels = 10,
                          pos_group_name = "Metastases", 
                          neg_group_name = "Encephalitis")

ggsave("./Figures/scRNA_Group_Volcano_plot_Metastases-Encephalitis.pdf", p_volcano, height = 6, width = 8)

# Identify group-associated differential genes relative to other groups.
Idents(scRNA) <- "Group"

CELLDEG <- FindAllMarkers(scRNA, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
CELLDEG <- CELLDEG[!is.na(CELLDEG$p_val), ]
write.csv(CELLDEG, "./Result/scRNA_Group_All_Markers.csv", row.names = FALSE)

# Export pseudobulk lollipop plots, group-specific CSVs and bar/dot plots.
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "00_Shared/plot_top_pseudobulk.R"))

plot_top_pseudobulk(scRNA, CELLDEG, group_by = "Group",sample_col = "HLA_Splitter",
                    top_n_lollipop = 20, top_n_bardot = 10, my_colors = GROUP_COLORS)

# Save analysis outputs.
saveRDS(scRNA, file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "CSF_scRNA.rds"))
