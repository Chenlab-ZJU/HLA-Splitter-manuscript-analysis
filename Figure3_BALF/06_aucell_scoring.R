# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 05_AUCell_scoring.R
# Data cohort: balf. Inputs use the compact English data layout.

# ==============================================================================
# File: 05_aucell_scoring.R
# Description: derive differential-gene signatures, score cells with AUCell and evaluate sample-level disease scores.
# ==============================================================================
rm(list = ls())
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure3_BALF/01_setup_and_helpers.R"))
# Read the saved object or source table used by this analysis.
Tcell <- readRDS(file.path("D:/项目/2026.Nature Communications/Revison/Data", "balf", "lavage_Tcell.rds"))
Sample_Info <- read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "balf", "SampleInfo.csv"))

# Retain target groups and drop unused factor levels.
cell_subset <- subset(Tcell, subset = Group %in% c("LC", "ILD", "CAP"))
cell_subset$HLA_Splitter %>% table()
cell_subset$Group <- factor(cell_subset$Group, levels = c("CAP", "LC", "ILD"))
comparison <- list(c("CAP","LC"), c("CAP","ILD"), c("ILD","LC"))

# ==============================================================================
# Identify disease-associated genes for data-derived signatures.
# ==============================================================================
Idents(cell_subset) <- "Group"
message("Calculating differential genes for each group...")

# Construct the gene-set list required by AUCell.
#geneSets <- read.csv("./Result/MSigDB/Lung_disease_DEGs.csv")
#geneSets = as.list(geneSets)
#geneSets = lapply(geneSets,function(x)x[nchar(x)>0])

# geneSets <- split(filtered_markers$gene, filtered_markers$cluster)

Group_DEG = read.csv('./Result/Tcell_Group_DEG.csv')
filtered_markers = Group_DEG %>%
  filter(PB_Log2FC>0.1)  %>%
  group_by(Group_Col)
geneSets <- split(filtered_markers$gene, filtered_markers$Group_Col)


# ==============================================================================
# Calculate AUCell disease-signature scores.
# ==============================================================================
message("Building AUCell rankings and calculating scores...")
names(geneSets) <- paste0(names(geneSets), "_Signature")
counts <- GetAssayData(cell_subset, layer = "data")
cells_rankings <- AUCell_buildRankings(counts, plotStats = FALSE, splitByBlocks = TRUE)

# Calculate AUC using the specified aucMaxRank; typical ranks cover 5-10% of genes.
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.1)

# Store scores in the Seurat cell metadata.
signature_names <- names(geneSets)
for (sig in signature_names) {
  cell_subset[[sig]] <- as.numeric(getAUC(cells_AUC)[sig, ])
}

# ==============================================================================
# Compare disease-associated scores across cell subsets.
# ==============================================================================
# Plot signature scores on UMAP.
p_umap_sig <- FeaturePlot(cell_subset, features = signature_names, ncol = 3, pt.size = 0.5, reduction = "umap") & 
  scale_color_gradientn(colors = c("lightgrey", "firebrick3")) & theme_void()
ggsave("./Figures/Tcell_Signature_UMAP.pdf", p_umap_sig, height = 4, width = 12)

# Compare signature scores across cell subtypes using violins.
p_vln_sig <- VlnPlot(cell_subset, group.by = "Celltype", features = signature_names, stack = TRUE, pt.size = 0, flip = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"), legend.position = "none")
ggsave("./Figures/Tcell_Signature_Celltype_Vln.pdf", p_vln_sig, height = 6, width = 8)

# Summarize cell-subtype disease scores in a dot plot.
p_dot_sig <- DotPlot(cell_subset, group.by = "Celltype", features = signature_names) +
  coord_flip() +
  scale_color_gradientn(colors = c("navy", "white", "firebrick3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(x = "Disease Signature", y = "Cell Type")
ggsave("./Figures/Tcell_Signature_Celltype_Dotplot.pdf", p_dot_sig, height = 4, width = 7)


# ==============================================================================
# Average cell scores by sample and evaluate ROC performance.
# ==============================================================================
# Extract sample-average scores.
df_sample <- FetchData(cell_subset, c("Group", "HLA_Splitter", signature_names))
df_agg <- aggregate(df_sample[, signature_names], by = list(Group = df_sample$Group, SampleName = df_sample$HLA_Splitter), FUN = mean)

# Draw sample-level score box plots.
df_melt <- reshape2::melt(df_agg, id.vars = c("Group", "SampleName"), variable.name = "Signature", value.name = "Score")
p_box_samp <- ggplot(df_melt, aes(x = Group, y = Score, fill = Group)) +
  facet_wrap(~ Signature, scales = "free_y") +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  stat_compare_means(comparisons = comparison, label = "p.signif") +
  my_theme()+
  scale_fill_manual(values = GLOBAL_COLORS) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(title = "Sample-level Disease Signature Scores", x = "Actual Group", y = "Mean AUCell Score")
p_box_samp
ggsave("./Figures/Tcell_Signature_Sample_Boxplot.pdf", p_box_samp, height = 5, width = 6)

# scatter plot (Sample)
df_agg_scale = df_agg
df_agg[,3:5] = apply(df_agg[,3:5],2,scale)
p_scatter_samp = ggplot(df_agg,aes(CAP_Signature, LC_Signature,color = Group,label = SampleName))+
  geom_point(size = 3, alpha = 0.5)+
  geom_abline(intercept = 0, slope = 1)+
  geom_text_repel(size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE,max.overlaps =30)+
  scale_color_manual(values = GLOBAL_COLORS) +
  my_theme()+ coord_fixed()+
  xlim(c(-2,3))+ylim(c(-2,3))
p_scatter_samp  
ggsave(paste0("./Figures/TCell_signature_scatter.pdf"),p_scatter_samp,height = 4,width = 5.5)

# Evaluate score discrimination using ROC curves.
Roc_All <- data.frame()
AUC_labels <- c()

# Calculate ROC curves for each of the three disease groups.
for (grp in levels(cell_subset$Group)) {
  sig_name <- paste0(grp, "_Signature")
  # Binary outcome: membership of the current disease group.
  true_labels <- ifelse(df_agg$Group == grp, 1, 0)
  scores <- df_agg[[sig_name]]
  
  # Calculate ROC only when the required group observations exist.
  if (sum(true_labels) > 0) {
    Roc <- plot.roc(true_labels, scores, percent = TRUE, quiet = TRUE)
    AUC_labels <- c(AUC_labels, paste0(grp, " Model AUC: ", round(Roc$auc, 2), "%"))
    
    tmp_roc_df <- data.frame(
      sensitivity = Roc$sensitivities,
      specificity = Roc$specificities,
      Model = sig_name
    )
    Roc_All <- rbind(Roc_All, tmp_roc_df)
  }
}
Roc_All <- Roc_All %>% 
  arrange(Model, sensitivity)

# Draw the group-specific ROC curves together.
p_roc <- ggplot(Roc_All, aes(x = specificity, y = sensitivity, color = Model)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 100, linetype = "dashed", color = "grey50") +
  theme_classic(base_size = 14) +
  scale_x_reverse() +
  labs(title = "Diagnostic Performance of AUCell Signatures",
       x = "Specificity (%)", y = "Sensitivity (%)") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = c(0.5, 0.25), 
        legend.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Add AUC text labels.
y_pos <- c(35, 30, 25)
for(i in seq_along(AUC_labels)) {
  p_roc <- p_roc + ggplot2::annotate("text", x = 30, y = y_pos[i], label = AUC_labels[i], size = 4.5, fontface = "bold")
}
p_roc
ggsave("./Figures/Tcell_Signature_ROC_Plot.pdf", p_roc, dpi = 300, height = 5, width = 5)

# Save the Seurat object containing the scores.
# Export the analysis source data for figure reproduction.
saveRDS(cell_subset, "./Result/lavage_Tcell_Scored.rds")

# Save disease-signature scores for each sample.
write.csv(df_agg, "./Result/Disease_signature_score.csv")
