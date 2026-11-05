# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream T-cell annotation, repertoire, differential-expression and balanced ligand-receptor analysis.
# Original source: 07_Cell_communication.R
# Data cohort: urine. Inputs use the compact English data layout.

# ==============================================================================
# File: 07_Cell_communication_SampleLevel.R
# Description: sample-level cell-communication analysis from single-cell data.
# Use HLA_Splitter as the sample identifier and compare communication across Group.
# ==============================================================================

# Confirm that the helper-script path exists.
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure4_urine_immune/01_setup_and_helpers.R"))
# ==============================================================================
# Configure parameters and initialize the environment.
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(tibble)
  library(circlize)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(grid) # Draw diagonal missing-value markers in ComplexHeatmap.
})

outdir <- "./Result/Cellchart_SampleLevel"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Split the analysis by the configured sample_col.
sample_col <- "HLA_Splitter"
group_col <- "Group"
celltype_col <- "Major_celltype"

target_group_order <- c("NR", "TCMR", "ABMR", "Mix", "BK")
comparison = list(c("NR","TCMR"),c("NR","ABMR"),c("NR","Mix"),c("NR","BK"))

organism <- "human"      
assay_use <- "RNA"       
min_cells <- 10          
n_workers <- 6

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 20 * 1024^3)
future::plan("multisession", workers = n_workers)

# ==============================================================================
# Define core analysis helpers.
# ==============================================================================

# Retrieve expression data compatibly with Seurat v4/v5.
get_data_matrix <- function(seu, assay = "RNA") {
  mat <- tryCatch({ 
    GetAssayData(seu, assay = assay, layer = "data") 
  }, error = function(e) { 
    GetAssayData(seu, assay = assay, slot = "data") 
  })
  return(mat)
}

# Run CellChat for one sample; the input must already be subset by sample.
run_cellchat_one_sample <- function(seu, sample_name, celltype_col, assay_use, organism, min_cells) {
  message("\n>> Constructing sample-level CellChat object: ", sample_name)
  
  DefaultAssay(seu) <- assay_use
  meta <- seu@meta.data
  
  celltypes_char <- as.character(meta[[celltype_col]])
  ct_n <- table(celltypes_char)
  keep_ct <- names(ct_n)[ct_n >= min_cells]
  
  if (length(keep_ct) < 2) {
    warning("Sample [", sample_name, "] has fewer than two eligible cell types; skipping.")
    return(NULL)
  }
  
  seu_sub <- subset(seu, cells = rownames(meta)[celltypes_char %in% keep_ct])
  data.input <- get_data_matrix(seu_sub, assay = assay_use)
  
  meta_sub <- seu_sub@meta.data
  meta_sub$samples <- factor(sample_name)
  meta_sub$celltype_for_CellChat <- factor(as.character(meta_sub[[celltype_col]]))
  
  # Align expression-matrix columns with metadata rows.
  common_cells <- intersect(colnames(data.input), rownames(meta_sub))
  data.input <- data.input[, common_cells, drop = FALSE]
  meta_sub <- meta_sub[common_cells, , drop = FALSE]
  
  if (ncol(data.input) == 0) {
    warning("Sample [", sample_name, "] has no cells available in the expression matrix; skipping.")
    return(NULL)
  }
  
  if (length(unique(meta_sub$celltype_for_CellChat)) < 2) {
    warning("Sample [", sample_name, "] has fewer than two eligible cell types after filtering; skipping.")
    return(NULL)
  }
  
  cellchat <- createCellChat(
    object = data.input,
    meta = meta_sub,
    group.by = "celltype_for_CellChat"
  )
  
  cellchat@DB <- if (organism == "human") {
    CellChatDB.human
  } else {
    CellChatDB.mouse
  }
  
  cellchat <- subsetData(cellchat)
  if (is.null(cellchat@data.signaling) || nrow(cellchat@data.signaling) == 0) {
    warning("Sample [", sample_name, "] has no available CellChat signaling genes; check gene-symbol identifiers.")
    return(NULL)
  }
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  if (is.null(cellchat@LR$LRsig) || nrow(cellchat@LR$LRsig) == 0) {
    warning("Sample [", sample_name, "] has no significant ligand-receptor pairs; skipping.")
    return(NULL)
  }
  
  cellchat <- tryCatch({
    computeCommunProb(
      cellchat,
      type = "triMean",
      raw.use = TRUE,
      population.size = TRUE
    )
  }, error = function(e) {
    warning("Sample [", sample_name, "] computeCommunProb failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(cellchat)) return(NULL)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  
  if (is.null(cellchat@net$prob) || length(cellchat@net$prob) == 0) {
    warning("Sample [", sample_name, "] has no interactions after filterCommunication; skipping.")
    return(NULL)
  }
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  if (!is.null(cellchat@netP$prob) && length(cellchat@netP$prob) > 0) {
    cellchat <- tryCatch({
      suppressMessages(
        netAnalysis_computeCentrality(cellchat, slot.name = "netP")
      )
    }, error = function(e) {
      message("Sample [", sample_name, "] centrality calculation failed; retaining the CellChat object: ", e$message)
      cellchat
    })
  }
  message(">> Completed sample-level CellChat object: ", sample_name)
  return(cellchat)
}

# Extract global interaction metrics.
extract_overall_interaction <- function(cellchat, sample_name, group_name) {
  count_mat <- cellchat@net$count
  weight_mat <- cellchat@net$weight
  data.frame(
    Sample = sample_name,
    Group = group_name,
    total_interaction_number = sum(count_mat),
    total_interaction_strength = sum(weight_mat)
  )
}

# Extract pathway information-flow metrics.
extract_pathway_information_flow <- function(cellchat, sample_name, group_name) {
  prob <- cellchat@netP$prob
  if (is.null(prob) || length(dim(prob)) < 3) return(data.frame())
  pathway_flow <- apply(prob, 3, sum)
  data.frame(
    Sample = sample_name,
    Group = group_name,
    pathway = names(pathway_flow),
    information_flow = as.numeric(pathway_flow)
  )
}

# Calculate between-group Wilcoxon p values from sample-level observations.
calculate_group_stats <- function(df_long, value_col, feature_col) {
  grps <- sort(unique(df_long$Group))
  res <- list()
  
  for (cmb in combn(grps, 2, simplify = FALSE)) {
    g1 <- cmb[1]; g2 <- cmb[2]
    
    # Complete sample-pathway combinations; absent detected pathways are assigned zero by this workflow.
    df_wide <- df_long %>%
      filter(Group %in% c(g1, g2)) %>%
      select(Sample, Group, !!sym(feature_col), !!sym(value_col)) %>%
      pivot_wider(names_from = !!sym(feature_col), values_from = !!sym(value_col), values_fill = 0) %>%
      pivot_longer(cols = -c(Sample, Group), names_to = feature_col, values_to = value_col)
    
    # Calculate summary statistics.
    stats_df <- df_wide %>%
      group_by(!!sym(feature_col)) %>%
      summarise(
        mean_g1 = mean(.data[[value_col]][Group == g1]),
        mean_g2 = mean(.data[[value_col]][Group == g2]),
        p_val = tryCatch({ wilcox.test(.data[[value_col]] ~ Group)$p.value }, error=function(e) NA),
        .groups = "drop"
      ) %>%
      mutate(
        comparison = paste0(g2, "_vs_", g1),
        log2FC = log2((mean_g2 + 1e-6) / (mean_g1 + 1e-6)),
        p_val_adj = p.adjust(p_val, method = "BH")
      ) %>%
      filter(!is.na(p_val)) # Exclude pathways that are zero in both groups and cannot be tested.
    
    res[[paste0(g2, "_vs_", g1)]] <- stats_df
  }
  bind_rows(res) %>% arrange(comparison, p_val)
}


# ==============================================================================
# Load data and construct metadata mappings.
# ==============================================================================
# Read the saved object or source table used by this analysis.
scRNA <- readRDS(find_object_file("Urine_scRNA.rds", "Urine_ALL.rds"))

# Clean the input data.
scRNA@meta.data[[group_col]] <- as.character(scRNA@meta.data[[group_col]])
scRNA@meta.data[[sample_col]] <- as.character(scRNA@meta.data[[sample_col]])

# Map each sample to its clinical group.
sample_meta <- scRNA@meta.data %>% 
  dplyr::select(all_of(sample_col), all_of(group_col)) %>% 
  distinct() %>% 
  drop_na()

sample_to_group <- setNames(sample_meta[[group_col]], sample_meta[[sample_col]])
sample_levels <- names(sample_to_group)

message("Found ", length(sample_levels), " independent samples across ", length(unique(sample_to_group)), " groups.")

# ==============================================================================
# Construct and save sample-level CellChat objects.
# ==============================================================================
cellchat_list <- list()
for (smp in sample_levels) {
  cells_to_keep <- rownames(scRNA@meta.data)[scRNA@meta.data[[sample_col]] == smp]
  
  if (length(cells_to_keep) == 0) {
    message("Warning: sample ", smp, " has no matching cells in scRNA.")
    next
  }
  
  seu_sub <- scRNA[, cells_to_keep]
  if (!sample_col %in% colnames(seu_sub@meta.data) | 
      !celltype_col %in% colnames(seu_sub@meta.data)) {
    message("Warning: sample ", smp, " is missing required metadata columns.")
    next
  }
  
  cellchat_one <- tryCatch({
    run_cellchat_one_sample(
      seu = seu_sub,
      sample_name = smp,
      celltype_col = celltype_col,
      assay_use = assay_use,
      organism = organism,
      min_cells = min_cells
    )
  }, error = function(e) {
    message("Error: sample ", smp, " CellChat analysis failed: ", e$message)
    return(NULL)
  })
  
  if (!is.null(cellchat_one)) {
    cellchat_list[[smp]] <- cellchat_one
  }
  gc()
}

cellchat_list <- cellchat_list[!sapply(cellchat_list, is.null)]

if (length(cellchat_list) < 2) {
  stop("Fewer than two eligible CellChat samples remain; group comparisons cannot proceed.")
}

# Export the analysis source data for figure reproduction.
saveRDS(cellchat_list, file.path(outdir, "01_CellChat_Sample_objects.rds"))


# ==============================================================================
# Extract summaries and test group differences.
# ==============================================================================

# Plot overall communication summaries as box plots.
overall_summary <- bind_rows(lapply(names(cellchat_list), function(s) {
  extract_overall_interaction(cellchat_list[[s]], sample_name = s, group_name = sample_to_group[s])
}))
write.csv(overall_summary, file.path(outdir, "02_SampleLevel_Overall_Stats.csv"), row.names = FALSE)

# Compare pathway information flow across groups.
pathway_flow_long <- bind_rows(lapply(names(cellchat_list), function(s) {
  extract_pathway_information_flow(cellchat_list[[s]], sample_name = s, group_name = sample_to_group[s])
}))
write.csv(pathway_flow_long, file.path(outdir, "03_SampleLevel_PathwayFlow_long.csv"), row.names = FALSE)

# Run the configured pairwise statistical comparisons.
pathway_diff_stats <- calculate_group_stats(pathway_flow_long, value_col = "information_flow", feature_col = "pathway")
write.csv(pathway_diff_stats, file.path(outdir, "04_GroupLevel_Pathway_Diff_Stats.csv"), row.names = FALSE)


# ==============================================================================
# Draw comparison figures with a consistent style.
# ==============================================================================

# Compare sample-level global communication strengths and significance.
actual_groups <- intersect(target_group_order, unique(overall_summary$Group))
overall_summary$Group = factor(overall_summary$Group, levels = actual_groups)

p_num <- ggboxplot(
  overall_summary, 
  x = "Group", 
  y = "total_interaction_number", 
  color = "Group",             # Match borders to fill colors rather than using black outlines.
  fill = "Group", alpha = 0.2, # Use opacity 0.2 for the highlight layer.
  palette = GROUP_COLORS,      # Apply the global named mapping.
  add = "jitter", 
  add.params = list(size = 1.5, alpha = 0.8, jitter = 0.2)
) +
  theme_classic(base_size = 14) +
  stat_compare_means(
    comparisons = comparison,  # Apply the global statistical-test configuration.
    method = "wilcox.test", 
    label = "p.format", 
    step.increase = 0.12,
    tip.length = 0.02,         # Refine line joins.
    bracket.size = 0.3         # Adjust comparison-bracket lines.
  ) +
  labs(x = NULL, y = "Interaction Counts") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )

p_str <- ggboxplot(
  overall_summary, 
  x = "Group", 
  y = "total_interaction_strength", 
  color = "Group",
  fill = "Group", alpha = 0.2,
  palette = GROUP_COLORS,      
  add = "jitter", 
  add.params = list(size = 1.5, alpha = 0.8, jitter = 0.2)
) +
  theme_classic(base_size = 14) +
  stat_compare_means(
    comparisons = comparison,    
    method = "wilcox.test", 
    label = "p.format",
    step.increase = 0.12,
    tip.length = 0.02,
    bracket.size = 0.3
  ) +
  labs(x = NULL, y = "Interaction Strength") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )

ggsave(fig_file("CellChat/Overview", 1, "Total_Interaction_Number_Strength"), p_num+p_str, width = 8, height = 4)

# Draw group box plots for top differential pathways.
top_pathway_dir <- ensure_dir(figure_path("CellChat/Top_Differential_Pathways_Boxplots"))
dir.create(top_pathway_dir, recursive = TRUE, showWarnings = FALSE)

for (cmp in unique(pathway_diff_stats$comparison)) {
  # Select differential pathways and complete missing observations with zeros.
  top_pwys <- pathway_diff_stats %>%
    filter(comparison == cmp & p_val < 1) %>%
    arrange(desc(abs(log2FC))) %>%
    slice_head(n = 6) %>%
    pull(pathway)
  
  if (length(top_pwys) == 0) {
    message("Comparison [", cmp, "] has no significant differential pathways; skipping the plot.")
    next
  }
  
  g_pair <- strsplit(cmp, "_vs_")[[1]]
  
  plot_df <- pathway_flow_long %>%
    filter(Group %in% g_pair & pathway %in% top_pwys) %>%
    select(Sample, Group, pathway, information_flow) %>%
    pivot_wider(names_from = pathway, values_from = information_flow, values_fill = 0) %>%
    pivot_longer(cols = all_of(top_pwys), names_to = "pathway", values_to = "information_flow") %>%
    mutate(pathway = factor(pathway, levels = top_pwys))
  
  # Draw consistently styled faceted box plots.
  p_pwys <- ggboxplot(
    plot_df, x = "Group", y = "information_flow", 
    color = "Group", fill = "Group", alpha = 0.2, 
    palette = GROUP_COLORS, add = "jitter", 
    facet.by = "pathway", scales = "free_y",
    add.params = list(size = 1.5, alpha = 0.8, jitter = 0.2)
  ) +
    stat_compare_means(
      comparisons = comparison,
      method = "wilcox.test", 
      label = "p.format", 
      step.increase = 0.12,
      tip.length = 0.02,
      bracket.size = 0.3
    ) +
    theme_classic(base_size = 12) +
    labs(title = paste0("Signaling Information Flow: ", cmp), x = NULL, y = "Information Flow") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6),
      strip.background = element_rect(fill = "#F5F5F5", color = NA), # Style facet-strip backgrounds.
      strip.text = element_text(face = "bold", size = 11)
    )
  
  cmp_i <- match(cmp, unique(pathway_diff_stats$comparison))
  ggsave(file.path(top_pathway_dir, sprintf("%02d_Pathway_Boxplot_%s.pdf", cmp_i, cmp)), p_pwys, width = 10, height = 6)
}

# Visualize global pathway strengths as a heatmap.
heatmap_data_long <- pathway_flow_long %>%
  group_by(Group, pathway) %>%
  summarise(mean_flow = mean(information_flow), .groups = "drop")

heatmap_mat <- heatmap_data_long %>%
  pivot_wider(names_from = Group, values_from = mean_flow, values_fill = 0) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

group_order <- intersect(target_group_order, colnames(heatmap_mat))
heatmap_mat <- heatmap_mat[, group_order, drop = FALSE]

row_sums <- rowSums(heatmap_mat)
heatmap_mat <- heatmap_mat[row_sums > 1e-3, , drop = FALSE]

heatmap_pdf <- fig_file("CellChat/Pathway_Heatmaps", 1, "Average_Pathway_Strength")
pdf_height <- max(6, nrow(heatmap_mat) * 0.15)

pdf(heatmap_pdf, width = 6, height = pdf_height)
pheatmap(
  heatmap_mat,
  scale = "row",           
  cluster_cols = FALSE,    
  cluster_rows = TRUE,     
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Average Pathway Strength (Z-score)",
  fontsize_row = 8,
  angle_col = "45"
)
dev.off()

# ==============================================================================
# Compare selected pathways using grouped bars.
# ==============================================================================
#target_pathways <- c("MHC-I", "MHC-II", "CXCL", "CCL", "TGFb", "TNF")
target_pathways <- unique(pathway_flow_long$pathway)

plot_data_all <- pathway_flow_long %>%
  select(Sample, Group, pathway, information_flow) %>%
  pivot_wider(names_from = pathway, values_from = information_flow, values_fill = 0) %>%
  pivot_longer(cols = -c(Sample, Group), names_to = "pathway", values_to = "information_flow") %>%
  mutate(Group = factor(Group, levels = group_order))

bar_dir <- ensure_dir(figure_path("CellChat/Target_Pathways_Barplots"))
dir.create(bar_dir, recursive = TRUE, showWarnings = FALSE)

for (pwy in target_pathways) {
  df_pwy <- plot_data_all %>% filter(pathway == pwy)
  if(nrow(df_pwy) == 0 || max(df_pwy$information_flow) == 0) {
    message("  - Pathway [", pwy, "] has no expression across samples; skipping the plot.")
    next
  }
  
  p_bar <- ggbarplot(
    df_pwy,
    x = "Group", 
    y = "information_flow",
    color = "Group",             # Use mapped outline colors instead of black borders.
    fill = "Group", alpha = 0.2, # Use semitransparent highlight fills.
    palette = GROUP_COLORS,      
    width = 0.65,                # Use narrower bars for clearer spacing.
    add = c("mean_se", "jitter"),
    error.plot = "upper_errorbar", # Show upward error bars only to avoid crossing sample points.
    add.params = list(
      size = 1.5, alpha = 0.8, jitter = 0.2,        
      linewidth = 0.6 
    ) ) +
    stat_compare_means(
      comparisons = comparison,    
      method = "wilcox.test", 
      label = "p.signif",          
      step.increase = 0.12,
      tip.length = 0.02,           # Refine connecting lines.
      bracket.size = 0.3
    ) +
    theme_classic(base_size = 14) +
    labs(title = paste0(pwy), x = NULL, y = "Signaling Strength") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6)
    )
  
  pwy_i <- match(pwy, target_pathways)
  ggsave(file.path(bar_dir, sprintf("%02d_Pathway_Barplot_%s.pdf", pwy_i, pwy)), p_bar, width = 3, height = 4)
}

# ==============================================================================
# Align combined heatmaps and mark missing values with diagonal lines.
# ==============================================================================
message("\n>> Extracting and plotting group heatmaps with aligned axes...")

extract_cell_pathway_strength <- function(cellchat, sample_name, group_name) {
  if (is.null(cellchat@netP$centr)) return(NULL)
  pwys <- names(cellchat@netP$centr)
  if (length(pwys) == 0) return(NULL)
  
  res_list <- lapply(pwys, function(pwy) {
    out_deg <- cellchat@netP$centr[[pwy]]$outdeg
    in_deg <- cellchat@netP$centr[[pwy]]$indeg
    if(is.null(out_deg) || is.null(in_deg)) return(NULL)
    data.frame(Sample = sample_name, Group = group_name, Pathway = pwy,
               CellType = names(out_deg), Outgoing = as.numeric(out_deg),
               Incoming = as.numeric(in_deg), stringsAsFactors = FALSE)
  })
  return(bind_rows(res_list))
}

cell_pwy_df <- bind_rows(lapply(names(cellchat_list), function(s) {
  extract_cell_pathway_strength(cellchat_list[[s]], sample_name = s, group_name = sample_to_group[s])
}))

cell_pwy_mean <- cell_pwy_df %>%
  group_by(Group, Pathway, CellType) %>%
  summarise(mean_Outgoing = mean(Outgoing, na.rm = TRUE),
            mean_Incoming = mean(Incoming, na.rm = TRUE), .groups = "drop")

heatmap_dir <- ensure_dir(figure_path("CellChat/GroupLevel_Combined_Heatmaps"))
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

all_cells <- sort(unique(cell_pwy_mean$CellType))

global_pwy_stat <- cell_pwy_mean %>%
  group_by(Pathway) %>%
  summarise(
    global_total_out = sum(mean_Outgoing, na.rm = TRUE),
    global_total_in = sum(mean_Incoming, na.rm = TRUE),
    global_max_out = max(mean_Outgoing, na.rm = TRUE), 
    global_max_in = max(mean_Incoming, na.rm = TRUE)   
  ) %>%
  filter(global_total_out > 0.01 & global_total_in > 0.01) %>%
  arrange(desc(global_total_out)) 

all_pwys <- global_pwy_stat$Pathway
max_out_dict <- setNames(global_pwy_stat$global_max_out, global_pwy_stat$Pathway)
max_in_dict <- setNames(global_pwy_stat$global_max_in, global_pwy_stat$Pathway)

group_pwy_totals <- cell_pwy_mean %>%
  group_by(Group, Pathway) %>%
  summarise(grp_total = sum(mean_Outgoing, na.rm = TRUE), .groups = "drop")
max_bar_val <- max(group_pwy_totals$grp_total, na.rm = TRUE)

col_fun <- circlize::colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                                c("white", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"))

actual_groups <- intersect(target_group_order, unique(cell_pwy_mean$Group))

for (grp_i in seq_along(actual_groups)) {
  grp <- actual_groups[grp_i]
  message("  - Drawing ", grp, " group combined heatmap with aligned axes...")
  df_grp <- cell_pwy_mean %>% filter(Group == grp)
  
  # Use tidyr::complete() to insert missing rows as NA.
  df_out_wide <- df_grp %>%
    select(Pathway, CellType, mean_Outgoing) %>%
    filter(Pathway %in% all_pwys) %>%
    complete(Pathway = all_pwys, CellType = all_cells, fill = list(mean_Outgoing = NA)) %>%
    pivot_wider(names_from = CellType, values_from = mean_Outgoing) %>%
    column_to_rownames("Pathway")
  mat_out <- as.matrix(df_out_wide)[all_pwys, all_cells] 
  
  df_in_wide <- df_grp %>%
    select(Pathway, CellType, mean_Incoming) %>%
    filter(Pathway %in% all_pwys) %>%
    complete(Pathway = all_pwys, CellType = all_cells, fill = list(mean_Incoming = NA)) %>%
    pivot_wider(names_from = CellType, values_from = mean_Incoming) %>%
    column_to_rownames("Pathway")
  mat_in <- as.matrix(df_in_wide)[all_pwys, all_cells]
  
  # Normalize using the global maximum.
  mat_out_norm <- mat_out
  mat_in_norm <- mat_in
  
  for (p in all_pwys) {
    if (max_out_dict[p] > 0) mat_out_norm[p, ] <- mat_out[p, ] / max_out_dict[p]
    if (max_in_dict[p] > 0) mat_in_norm[p, ] <- mat_in[p, ] / max_in_dict[p]
  }
  
  pwy_strength <- rowSums(mat_out, na.rm = TRUE)
  
  # Draw the combined heatmaps.
  ha_left <- rowAnnotation(
    `Total\nStrength` = anno_barplot(
      pwy_strength,
      bar_width = 0.7,
      ylim = c(0, max_bar_val), 
      gp = gpar(fill = "darkgray", col = NA),
      width = unit(2.5, "cm"),
      axis_param = list(direction = "reverse")
    ),
    annotation_name_rot = 0
  )
  
  # Draw the sender heatmap and use cell_fun to mark NA values.
  ht_out <- Heatmap(
    mat_out_norm,
    name = "Sender", col = col_fun,
    na_col = "#F5F5F5",
    cluster_rows = FALSE, cluster_columns = FALSE, 
    show_row_names = FALSE,
    column_title = "Sender Cells (Outgoing)", column_names_rot = 45,
    row_title = "Signaling Pathways"
  )
  
  # Draw the receiver heatmap.
  ht_in <- Heatmap(
    mat_in_norm,
    name = "Receiver", col = col_fun,
    na_col = "#F5F5F5",
    cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = TRUE, row_names_side = "right",
    column_title = "Receiver Cells (Incoming)", column_names_rot = 45
  )
  
  ht_list <- ha_left + ht_out + ht_in
  
  pdf_height <- max(6, length(all_pwys) * 0.15 + 2)
  pdf_width <- max(10, length(all_cells) * 0.3 * 2 + 5)
  
  pdf(file.path(heatmap_dir, sprintf("%02d_Combined_Aligned_%s.pdf", grp_i, grp)), width = pdf_width, height = pdf_height)
  draw(ht_list, ht_gap = unit(0.5, "cm"))
  dev.off()
}
# ==============================================================================
# Compare key signaling pathways using cell-interaction chord diagrams.
# ==============================================================================
chord_dir <- ensure_dir(figure_path("CellChat/Chord_Diagrams"))
dir.create(chord_dir, recursive = TRUE, showWarnings = FALSE)

target_pathways <- c("MHC-I", "MHC-II", "CXCL", "CCL", "TGFb", "TNF")

# Define a consistent cell-type color mapping.
all_cells_global <- sort(unique(unlist(lapply(cellchat_list, function(x) levels(x@idents)))))
num_cells <- length(all_cells_global)

qual_cols <- brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector <- unlist(mapply(brewer.pal, qual_cols$maxcolors, rownames(qual_cols)))
cell_colors <- setNames(col_vector[1:num_cells], all_cells_global)

# Draw the selected pathway chord diagram for each group.
for (grp_i in seq_along(actual_groups)) {
  grp <- actual_groups[grp_i]
  samps <- names(sample_to_group)[sample_to_group == grp]
  valid_samps <- intersect(samps, names(cellchat_list))
  
  if (length(valid_samps) == 0) next
  
  for (pwy in target_pathways) {
    mat_list <- list()
    
    # Retrieve pathway interaction matrices for each sample in the group.
    for (s in valid_samps) {
      cc <- cellchat_list[[s]]
      if (!is.null(cc@netP$pathways) && pwy %in% cc@netP$pathways) {
        
        if (length(dim(cc@netP$prob)) == 3) {
          mat <- cc@netP$prob[, , pwy]
        } else {
          mat <- cc@netP$prob
        }
        
        template <- matrix(0, nrow = num_cells, ncol = num_cells, 
                           dimnames = list(all_cells_global, all_cells_global))
        
        # Insert each local matrix into the global cell-type template.
        cells_s <- rownames(mat)
        template[cells_s, cells_s] <- mat
        mat_list[[s]] <- template
      }
    }
    
    if (length(mat_list) == 0) {
      message("    (Group ", grp, " has no detected ", pwy, " pathway in any sample; skipping)")
      next
    }
    
    # Calculate the mean interaction matrix across samples in the group.
    mean_mat <- Reduce("+", mat_list) / length(valid_samps)
    
    # Exclude cell types with no interactions in the selected pathway.
    keep_idx <- rowSums(mean_mat) > 0 | colSums(mean_mat) > 0
    if (sum(keep_idx) < 2) next # Retain at least two participating cell-type nodes for drawing links.
    
    plot_mat <- mean_mat[keep_idx, keep_idx, drop = FALSE]
    
    message("  - Drawing [", grp, "] group [", pwy, "] pathway chord diagram...")
    
    pwy_i <- match(pwy, target_pathways)
    pdf_path <- file.path(chord_dir, sprintf("%02d_Chord_%s_%s.pdf", (grp_i - 1) * length(target_pathways) + pwy_i, grp, pwy))
    pdf(pdf_path, width = 7, height = 7)
    par(mar = c(1, 1, 3, 1)) # Preserve outer margins.
    
    # Draw the chord diagram with the interaction matrix passed explicitly as net.
    tryCatch({
      netVisual_chord_cell(
        net = plot_mat,
        color.use = cell_colors[rownames(plot_mat)]
      )
      title(paste0(pwy, " Signaling Network (", grp, ")"), line = 1, cex.main = 1.2, font.main = 2)
    }, error = function(e) {
      message("    Warning: plotting ", grp, " / ", pwy, " failed: ", e$message)
    })
    
    dev.off()
  }
}

future::plan("sequential")
