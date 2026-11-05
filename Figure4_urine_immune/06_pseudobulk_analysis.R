# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream urine processing and sample-level signature analysis.
# Original source: 09_Pseudobulk_analysis.R
# Data cohort: urine. Inputs use the compact English data layout.

# ==============================================================================
# File: 04_group_marker_analysis.R
# Description: configurable pseudobulk differential expression, enrichment and combined plots.
# ==============================================================================
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure4_urine_immune/01_setup_and_helpers.R"))
source_if_exists(file.path(HELPER_DIR, "plot_top_pseudobulk.R"), required = TRUE)

if (!exists("scRNA")) scRNA <- readRDS(find_object_file("Urine_scRNA.rds", "Urine_ALL.rds"))
if (!exists("Tcell")) Tcell <- readRDS(find_object_file("Urine_Tcell.rds"))
if (!exists("Myeloid")) Myeloid <- readRDS(find_object_file("Urine_Myeloid.rds"))

# ==============================================================================
# Configuration: update these parameters when changing the analysis object.
# ==============================================================================
current_obj   <- scRNA 
analysis_name <- "scRNA"     # Output directory and filename prefix, such as Tcell_Subset.

group_col     <- "Group"            # Clinical-group column in cell metadata.
sample_col    <- "HLA_Splitter"     # Sample-identifier column in cell metadata.

my_colors     <- GROUP_COLORS       # Color palette used in these plots.
my_comparisons <- list(c("NR","TCMR"), c("NR","ABMR"), c("NR","Mix")) # Statistical comparison pairs.
# ==============================================================================
# Export pseudobulk lollipop plots, group-specific CSVs and bar/dot plots.
# ==============================================================================
# Use distinct output paths to avoid overwriting other analyses.
out_dir <- result_path(analysis_name, "Group_Markers")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_top_pseudobulk(
  seurat_obj = current_obj, 
  group_by = group_col,  
  sample_col = sample_col,
  top_n_lollipop = 10, 
  top_n_bardot = 10, 
  my_colors = my_colors, 
  out_dir = out_dir
)

# Combine top differential genes across subsets.
Group_DEG <- data.frame()
unique_groups <- as.character(unique(current_obj@meta.data[[group_col]]))

for (grp in unique_groups) {
  csv_path <- file.path(out_dir, paste0("Marker_Stats_", grp, ".csv"))
  if (file.exists(csv_path)) { 
    top_gene_tmp <- read.csv(csv_path)
    Group_DEG <- rbind(Group_DEG, top_gene_tmp)
  }
}
# Export the analysis source data for figure reproduction.
write.csv(Group_DEG, file.path(out_dir, paste0(analysis_name, "_Top_DEGs.csv")), row.names = FALSE)

# ==============================================================================
# Draw the global pseudobulk heatmap.
# ==============================================================================
message("Generating global pseudobulk heatmap for: ", analysis_name)

# Select genes with the largest log2FC for each cell type.
top_heatmap_genes <- Group_DEG %>%
  dplyr::group_by(Group_Col) %>%
  dplyr::slice_max(order_by = PB_Log2FC, n = 5, with_ties = FALSE) %>%
  dplyr::pull(gene) %>%
  unique()


# Custom gene list for the global scRNA analysis.
top_heatmap_genes <- c(
  "GNLY", "NKG7", "MT1M", "MT1H", "MT1G",
  "GZMB", "CXCL9", "CXCL10", "ICAM1",
  "S100A8", "C1QA", "STAT1", "IFITM1", "HMOX1","IFNG","CD3E", "CD7",
  "STAT2", "IRF7", "PDCD1", "CXCL8", "IFIT2", "CXCL2", "CXCL3", "ICAM1"
)

# Custom gene list for T cells.
top_heatmap_genes <- c(
  "CD4", "CD8A","CXCR3", "CD27", "STAT1",
  "TNF", "IFNG", "IL21", "CCL4", 
  "GZMK", "GZMB", "NKG7",
  "HLA-A","HLA-B" ,"HLA-C","HLA-DRA","HLA-DPA","HLA-DQA",
  "HAVCR2","TIGIT",  "PDCD1", "CTLA4"
)

# Custom gene list for myeloid cells.
top_heatmap_genes <- c("C1QC", "NLRP3", "TREM2", "SPP1", "APOC1", "APOE",
                  "HLA-A","HLA-B" ,"HLA-C","HLA-DRA","HLA-DPA","HLA-DQA",
                  "CXCL9", "CXCL10", "IFNG", "TNF","CCL3",  "CCL5")


seu_avg <- AverageExpression(current_obj, group.by = group_col, assays = "RNA", return.seurat = TRUE)
seu_avg <- ScaleData(seu_avg, features = rownames(seu_avg))
gene_matrix_sc <- seu_avg@assays$RNA$scale.data %>% as.matrix()

genes_to_plot <- intersect(top_heatmap_genes, rownames(gene_matrix_sc))

library(pheatmap)
plot_heatmap_sc_bulk <- pheatmap(
  gene_matrix_sc[genes_to_plot, ], 
  cluster_rows = TRUE, 
  cluster_cols = FALSE, 
  angle_col = 45, 
  cellwidth = 40, 
  cellheight = 12 
)

ggsave(fig_file(file.path(analysis_name, "Pseudobulk"), 1, "Group_DEG_Heatmap"), plot_heatmap_sc_bulk, height = 12, width = 8)

# ==============================================================================
# Draw cell-type-specific violin plots and GO enrichment results.
# ==============================================================================
valid_groups <- as.character(unique(Group_DEG$Group_Col))

for (grp in valid_groups) {
  message("\n>> Analyzing subset: ", grp)
  
  # Draw violin plots.
  plot_top_violins(
    seu_obj = current_obj,             
    deg_data = Group_DEG,
    target_group = grp,
    out_dir = out_dir,
    comparisons_list = my_comparisons,   
    plot_colors = my_colors,
    group_col = "Group_Col",
    lfc_col = "PB_Log2FC",
    vln_group_by = group_col,
    lfc_threshold = 0.1,
    top_n = 8
  )
  
  # Run and plot GO enrichment.
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
# Combine representative-gene plots at single-cell and sample levels.
# ==============================================================================
# Representative genes for the global scRNA analysis.
target_genes <- c(
  # ABMR-associated endothelial injury, NK-cell and chemokine programs.
  "NKG7", "MT1H","CXCL9", "CXCL10",
  "GNLY", "GZMB", "ICAM1",
  # TCMR-associated macrophage/T-cell infiltration and tissue injury programs.
  "S100A8", "C1QA", "STAT1", "IFITM1", "HMOX1","IFNG","CD3E", "CD7",
  # BK-associated antiviral, exhaustion and chronic-inflammatory programs.
  "STAT2", "IRF7", "PDCD1", "CXCL8", "IFIT2", "CXCL2", "CXCL3", "ICAM1"
)

# Representative genes for T cells.
target_genes <- c(
  "CXCR3", "CD27", "STAT1", "TNF",
  "IL21", "HAVCR2", "PDCD1", "CTLA4"
)

# Representative genes for myeloid cells.
target_genes <- c("TREM2", "SPP1", "APOE",
                  "CXCL9", "CXCL10", "IFNG")

valid_genes <- intersect(target_genes, rownames(current_obj))
if (length(valid_genes) == 0) stop("No target genes were found in the expression matrix; check gene symbols and the input object.")

expr_data <- FetchData(current_obj, vars = c(valid_genes, group_col, sample_col))
expr_data <- expr_data[!is.na(expr_data[[group_col]]), ]

# ------------------------------------------------------------------------------
# Combine single-cell violin plots.
# ------------------------------------------------------------------------------
sc_data_long <- expr_data %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = valid_genes)) 

p_vln_combined <- plot_violin_comparison(seurat_obj = current_obj, 
                                         group_col = "Group", 
                                         comparison_groups = comparison,
                                         features_to_plot = valid_genes, 
                                         plot_colors = GLOBAL_COLORS)

ggsave(fig_file(file.path(analysis_name, "Pseudobulk"), 2, "Core_Genes_Single_Cell_Violin"), p_vln_combined, width = 12, height = 6)

# ------------------------------------------------------------------------------
# Combine sample-level mean-expression bar plots.
# ------------------------------------------------------------------------------
pb_data_long <- expr_data %>%
  group_by(!!sym(group_col), !!sym(sample_col)) %>%
  summarise(across(all_of(valid_genes), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = all_of(valid_genes), names_to = "Gene", values_to = "Mean_Expression") %>%
  mutate(Gene = factor(Gene, levels = valid_genes))

p_bar_combined <- ggplot(pb_data_long, aes(x = .data[[group_col]], y = Mean_Expression)) +
  stat_summary(
    fun = mean, geom = "bar", 
    aes(fill = .data[[group_col]], color = .data[[group_col]]),
    alpha = 0.5, width = 0.65, linewidth = 0.6
  ) +
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      data.frame(y = m, ymin = m, ymax = m + sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
    },
    geom = "errorbar", aes(color = .data[[group_col]]), width = 0.2, linewidth = 0.7
  ) +
  geom_jitter(aes(color = .data[[group_col]]), width = 0.18, size = 1.5, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_color_manual(values = my_colors, aesthetics = c("color", "fill")) +
  theme_classic(base_size = 14) +
  stat_compare_means(
    comparisons = my_comparisons, method = "wilcox.test", 
    label = "p.signif", step.increase = 0.12, tip.length = 0.02, bracket.size = 0.3
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
  labs(
    title = paste0(analysis_name, " - Pseudo-bulk Mean Expression"), 
    x = NULL, 
    y = "Mean Expression"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 15)),
    strip.text = element_text(face = "bold.italic", size = 14), 
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )

ggsave(fig_file(file.path(analysis_name, "Pseudobulk"), 3, "Core_Genes_Pseudobulk_Barplot"), p_bar_combined, width = 8, height = 12)
