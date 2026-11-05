# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis.
# Original source: 05_Tcell_subclustering.R
# Data cohort: urine. Inputs use the compact English data layout.

# ==============================================================================
# File: 04_celltype_subclustering.R
# Description: recluster myeloid/T-cell subsets, identify markers and score polarization programs.
# ==============================================================================
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure4_urine_immune/01_setup_and_helpers.R"))
#scRNA <- readRDS("./Result/Objects/Urine_scRNA.rds")
# Read the saved object or source table used by this analysis.
Tcell = readRDS(file.path("D:/项目/2026.Nature Communications/Revison/Data", "urine", "Urine_Tcell.rds"))

# T-cell analysis.
Tcell = scRNA[,scRNA$Major_celltype %in% c("T cell", "NK-like", "Cycling")]

Tcell <- Tcell %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 800) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, min.dist = 0.3) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8)

# Identify and visualize global marker genes.
export_top_markers(Tcell, "seurat_clusters", "Tcell")

# Manual checkpoint: fill Celltype in Result/Tcell/Markers/Tcell_makers_summary.csv before continuing.
df <- read.csv("./Result/Tcell/Markers/Tcell_makers_summary.csv")
Celltype <- df$Celltype
names(Celltype) <- as.character(df$seurat_clusters)
Tcell <- SetIdent(Tcell, value = "seurat_clusters")
Tcell <- RenameIdents(Tcell, Celltype)
Tcell$Celltype <- Idents(Tcell)
unique(Tcell$Celltype) %>% print()

# Exclude contaminant cell populations.
Tcell = Tcell[,!Tcell$Celltype %in% c("Contament", "T_stress", "T_MT")]
Tcell$Celltype <- factor(Tcell$Celltype, levels = c("NK","CD4_Tcm","T_Quiescent","CD4_Tph","Treg","T_Activated",
                                                    "CD8_Teff","CD8_Trm","T_ISG", "T_cycling"))
Tcell = RunUMAP(Tcell, reduction = "harmony", dims = 1:30, min.dist = 0.3)

# Draw the three-panel global UMAP.
p1 <- DimPlot(Tcell, group.by = "seurat_clusters", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p2 <- DimPlot(Tcell, group.by = "Celltype", label = TRUE, alpha = 0.5) + scale_color_manual(values = GLOBAL_COLORS)
p3 <- DimPlot(Tcell, group.by = "Group", alpha = 0.5) + scale_color_manual(values = GROUP_COLORS)
ggsave(fig_file("Tcell/UMAP_Markers", 1, "UMAP_Cluster_Celltype_Group"), p1 + p2 + p3, height = 5, width = 16)

p_F = FeaturePlot(Tcell, features = c("CXCR3", "CXCR4", "IFNG", "FCGR3A"),ncol=2,cols = c("grey","firebrick"))
ggsave(fig_file("Tcell/UMAP_Markers", 2, "FeaturePlot_CXCR3_CXCR4_IFNG_FCGR3A"), p_F, height = 8, width = 8)

FeaturePlot(Tcell, features = c("CD3D", "KLRD1"),ncol=2,cols = c("grey","firebrick"))
# Plot global marker genes.
global_mk <- c("CCR7", "LEF1", "CD27", "IGFL2",
               "KLRB1","IL21", "CD40LG", "FOXP3", "IL2RA",
               "CCR6", "IL23R", "GZMK", "GZMB" ,
               "CXCR4", "IFNG","IFIT3", "ISG15","MKI67", "TYMS","KLRD1","FCGR3A")
p_dot <- DotPlot(Tcell, group.by = "Celltype", features = global_mk) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + labs(x=NULL, y=NULL) +
  scale_color_gradientn(colours = c("navy","white","firebrick3"))
ggsave(fig_file("Tcell/UMAP_Markers", 3, "DotPlot_Celltype_Marker_Genes"), plot = p_dot, height = 4, width = 9)

# Draw group-faceted UMAPs and cell-proportion bars.
combined_plot <- ggarrange(plotlist = plot_grouped_umap(Tcell, "Celltype", "Group", GLOBAL_COLORS), 
                           common.legend = TRUE, legend = "right")
ggsave(fig_file("Tcell/UMAP_Markers", 4, "Grouped_UMAP_by_Rejection_Group"), combined_plot, height = 6, width = 8)

comparison <- list(c("NR","TCMR"),c("NR","ABMR"),c("NR","Mix"))
p_list = plot_celltype_ratio(seurat_obj = Tcell,
                             sample_col = "Sample.ID",group_col = "Group",
                             cluster_col = "Celltype",
                             comparison_groups = comparison,
                             group_colors = GROUP_COLORS, 
                             cluster_colors = GLOBAL_COLORS)

ggsave(fig_file("Tcell/Celltype_Composition", 1, "Celltype_Stacked_Bar_by_Sample"), p_list$stacked_bar, width = 8, height = 4)
ggsave(fig_file("Tcell/Celltype_Composition", 2, "Celltype_Proportion_Heatmap"), p_list$heatmap, width = 5, height = 5)

# Plot cell-subtype proportions.
comparison <- list(c("NR","TCMR"),c("NR","ABMR"),c("NR","Mix"))
p_cell_ratio <- plot_cell_ratio_comparison(cell_subset = Tcell,
                                           sample_col = "HLA_Splitter", cluster_col =  "Celltype", group_col = "Group",
                                           facet = T, comparison_groups = comparison, 
                                           write_to_csv = "./Result/Tcell_celltype_ratio", plot_colors = GLOBAL_COLORS)
ggsave(fig_file("Tcell/Celltype_Composition", 3, "Celltype_Proportion_Group_Comparison"), p_cell_ratio, width = 10, height = 10)
# Export the analysis source data for figure reproduction.
saveRDS(Tcell, file.path("D:/项目/2026.Nature Communications/Revison/Data", "urine", "Urine_Tcell.rds"))


# ==============================================================================
# Identify disease-group differential genes in T cells.
# ==============================================================================
Idents(Tcell) <- "Group"
ident1 = "TCMR"
ident2 = "NR"
# Discover marker genes with the original assay and comparison settings.
CELLDEG <- FindMarkers(Tcell, ident.1 =ident1, ident.2 = ident2, logfc.threshold = 0.1, min.pct = 0.1)

noise_pattern <- "^MT-|^HB[ABDEGQZ]|^RP[0-9]+-|^AC[0-9]+\\.|^AL[0-9]+\\.|^AP[0-9]+\\.|^CT[ABCD]-|^RPS|^RPL|MALAT1|NEAT1"
CELLDEG$gene = rownames(CELLDEG)
CELLDEG <- CELLDEG[!is.na(CELLDEG$p_val),] %>%
filter(!grepl(noise_pattern, gene))

CELLDEG$gene <- rownames(CELLDEG)
write.csv(CELLDEG, paste0("./Result/Tcell/Differential_Expression/Tcell_Group_marker_",ident1,"-",ident2,".csv"))

# Plot differential-expression volcanoes and GO enrichment.
p_volcano <- plot_volcano(diffexp = CELLDEG,
                          top_n_labels = 10,
                          pos_group_name = ident1, 
                          neg_group_name = ident2)

ggsave(fig_file("Tcell/Differential_Expression", 1, paste0("Volcano_", ident1, "_vs_", ident2)), p_volcano, height = 5, width = 6)

# Compare gene expression.
selected_gene = c("CXCR3","CXCR4","CD27", "IFNG", "GZMK", "GZMB" ,"PDCD1","CTLA4")
selected_gene = c("CXCR3", "CD27",
                  "TNF", "IFNG",
                  "HLA-A","HLA-B" ,"HLA-C","HLA-DRA","HLA-DPA","HLA-DQA",
                  "HAVCR2","TIGIT",  "PDCD1", "CTLA4")

p1 = plot_violin_comparison(seurat_obj = Tcell, group_col = "Group", 
                            comparison_groups = comparison, 
                            features_to_plot = selected_gene, plot_colors = GLOBAL_COLORS)
ggsave(fig_file("Tcell/Differential_Expression", 2, "Marker_Genes_Group_Violin"), p1, height = 9, width = 8)

# Save analysis outputs.
saveRDS(Tcell,file.path("D:/项目/2026.Nature Communications/Revison/Data", "urine", "Urine_Tcell.rds"))


# ==============================================================================
# Examine CXCR3 expression enrichment.
# ==============================================================================
target_gene <- "CXCR3"
TCMR_subset <- subset(Tcell, subset = Group == "TCMR")

p_Vln_TCMR <- VlnPlot(TCMR_subset, features = target_gene, group.by = "Celltype", pt.size = 0) + 
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, linewidth = 0.6) + # Show the median explicitly in box plots.
  theme_classic(base_size = 14) +
  labs(
    title = bquote(italic(.(target_gene))~"Expression in TCMR Group"), 
    x = NULL, 
    y = "Expression Level"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "none"
  )

# Calculate the CXCR3-positive fraction for each cell type within TCMR.
expr_data_tcmr <- FetchData(TCMR_subset, vars = c(target_gene, "Celltype"), layer = "data")
colnames(expr_data_tcmr) <- c("Expression", "Cluster")

cxcr3_tcmr_stats <- expr_data_tcmr %>%
  group_by(Cluster) %>%
  summarise(
    Total_Cells = n(),                                     
    Pos_Cells = sum(Expression > 0),                       
    Mean_Expression = mean(Expression),                    
    Mean_Pos_Expression = mean(Expression[Expression > 0]),
    Pos_Ratio = (Pos_Cells / Total_Cells) * 100            
  ) %>%
  arrange(desc(Pos_Ratio)) # Order cell types by descending CXCR3-positive fraction.

print("--- CXCR3-positive fractions by cell type within TCMR ---")
print(cxcr3_tcmr_stats)

p_bar_TCMR <- ggplot(cxcr3_tcmr_stats, aes(x = reorder(Cluster, -Pos_Ratio), y = Pos_Ratio, fill = Cluster)) +
  geom_col(color = "black", linewidth = 0.6, width = 0.7) + 
  geom_text(aes(label = sprintf("%.1f%%", Pos_Ratio)), vjust = -0.8, size = 3.5, fontface = "bold") +
  theme_classic(base_size = 14) +
  labs(
    title = bquote("Proportion of "~italic(.(target_gene))^"+"~" Cells in TCMR"), 
    x = "Cell Type in TCMR",
    y = "Percentage of Positive Cells (%)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

final_plot <- p_Vln_TCMR + p_bar_TCMR
print(final_plot)

ggsave(fig_file("Tcell/CXCR3", 1, "TCMR_CXCR3_Expression_by_Celltype"), final_plot, height = 9, width = 8)

# ==============================================================================
# Compare cytotoxic-gene expression between TCMR and ABMR.
# ==============================================================================
cytox_genes <- c("GZMA", "GZMB", "GZMK", "PRF1", "GNLY", "NKG7", "FASLG")
valid_cytox_genes <- intersect(cytox_genes, rownames(Tcell))

# Retain TCMR and ABMR cells only.
rejection_subset <- subset(Tcell, subset = Group %in% c("TCMR", "ABMR"))
Idents(rejection_subset) <- "Group"

# Test group differences for every specified cytotoxic gene.
cytox_deg <- FindMarkers(
  rejection_subset, 
  ident.1 = "TCMR", 
  ident.2 = "ABMR", 
  features = valid_cytox_genes,
  logfc.threshold = 0,
  min.pct = 0
) %>%
  mutate(Gene = rownames(.)) %>%
  relocate(Gene) %>%
  arrange(p_val_adj)

print("--- Cytotoxic-gene differential expression: positive values indicate higher expression in TCMR than ABMR ---")
print(cytox_deg)

# Plot core cytotoxic genes by clinical group using a dot plot.
p_cytox_dot <- DotPlot(rejection_subset, features = valid_cytox_genes, group.by = "Group") +
  theme_classic(base_size = 14) +
  labs(
    title = "Cytotoxic Gene Expression: TCMR vs ABMR",
    x = "Cytotoxic Effectors", 
    y = "Rejection Type"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic", color = "black"), # Italicize gene symbols.
    axis.text.y = element_text(color = "black", face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.line = element_blank()
  ) +
  scale_color_gradientn(colours = colorRampPalette(c("navy", "white", "firebrick3"))(100))

# Draw stacked split violin plots.
stacked_vln_plot = plot_violin_comparison(seurat_obj = rejection_subset, 
                                          group_col = "Group", 
                                          comparison_groups = comparison, 
                                          features_to_plot = c("GZMB", "GZMK", "PRF1","NKG7"), 
                                          plot_colors = GLOBAL_COLORS,
                                          ncol = 4)

# Combine and export the dot plot and stacked violins.
combine_plot <- cowplot::plot_grid(
  p_cytox_dot, 
  stacked_vln_plot, 
  ncol = 1, 
  labels = c("A", "B"),
  rel_heights = c(2, 2)  
)

print(combine_plot)
ggsave(fig_file("Tcell/Differential_Expression", 3, "TCMR_vs_ABMR_Cytotoxic_Genes"), combine_plot, height = 6, width = 8)
