# Purpose: Transitive shared function dependency.
# Original source: plot_volcano.R
# Data cohort: shared. Inputs use the compact English data layout.

#' Plot a highly optimized volcano plot for two-group comparison
#'
#' @param diffexp A data.frame of differential expression results containing: 'p_val_adj', 'avg_log2FC', and 'gene'.
#' @param pos_group_name Character, name of the group with positive log2FC (e.g., "Treatment").
#' @param neg_group_name Character, name of the group with negative log2FC (e.g., "Control").
#' @param p_val_adj_threshold The adjusted p-value threshold (default: 0.05).
#' @param log2fc_threshold The log2 fold change threshold (default: 1).
#' @param colors A named character vector for coloring points.
#' @param top_n_labels Number of top genes to label per direction (default: 10, means 10 Up and 10 Down).
#' @param title Title of the plot.
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_size_manual scale_alpha_manual geom_vline geom_hline labs theme element_text element_blank element_rect ggtitle
#' @importFrom ggrepel geom_text_repel

plot_volcano <- function(diffexp,
                        pos_group_name = "Group A",
                        neg_group_name = "Group B",
                        p_val_adj_threshold = 0.05,
                        log2fc_threshold = 1,
                        colors = c('Down' = "#3498db", 'Un-sigif' = "grey80", 'Up' = "#e74c3c"),
                        top_n_labels = 10,
                        title = "Differential Expression") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required.")
  
  # Replace zero adjusted p values with the smallest positive floating-point value.
  diffexp$p_val_adj_capped <- ifelse(diffexp$p_val_adj == 0, .Machine$double.xmin, diffexp$p_val_adj)
  diffexp$logP <- -log10(diffexp$p_val_adj_capped)
  
  # Label upregulated and downregulated genes.
  diffexp$change <- "Un-sigif"
  diffexp$change[diffexp$p_val_adj < p_val_adj_threshold & diffexp$avg_log2FC >= log2fc_threshold] <- "Up"
  diffexp$change[diffexp$p_val_adj < p_val_adj_threshold & diffexp$avg_log2FC <= -log2fc_threshold] <- "Down"
  
  # Set factor levels so nonsignificant points are drawn below significant points.
  diffexp$change <- factor(diffexp$change, levels = c("Un-sigif", "Down", "Up"))
  
  # Select the top differential genes for labels.
  # Rank genes using expression prevalence, absolute log2FC and -log10(p).
  if(all(c("pct.1", "pct.2") %in% colnames(diffexp))) {
    diffexp$base_exp <- (diffexp$pct.1 + diffexp$pct.2) / 2
  } else {
    diffexp$base_exp <- 1 
  }
  
  # Constant ranking terms contribute zero instead of producing NaN scores.
  safe_zscore <- function(x) {
    result <- rep(0, length(x))
    valid <- is.finite(x)
    if (sum(valid) > 1L && stats::sd(x[valid]) > 0) {
      result[valid] <- as.numeric(scale(x[valid]))
    }
    result
  }
  # Calculate a standardized composite ranking score.
  diffexp$composite_score <- 0
  sig_idx <- diffexp$change %in% c("Up", "Down")
  
  if(sum(sig_idx) > 0) {
    diffexp$composite_score[sig_idx] <- 
      safe_zscore(-log10(diffexp$p_val_adj_capped[sig_idx])) +  # Significance contribution.
      safe_zscore(abs(diffexp$avg_log2FC[sig_idx])) +           # Effect-size contribution.
      safe_zscore(diffexp$base_exp[sig_idx])                    # Expression-prevalence contribution.
  }
  
  # Rank upregulated and downregulated genes separately by composite score.
  top_up_genes <- diffexp[diffexp$change == "Up", ]
  top_up_genes <- top_up_genes[order(top_up_genes$composite_score, decreasing = TRUE), ]
  
  top_down_genes <- diffexp[diffexp$change == "Down", ]
  top_down_genes <- top_down_genes[order(top_down_genes$composite_score, decreasing = TRUE), ]
  
  # Select the top N genes for labels.
  label_genes_df <- rbind(head(top_up_genes, top_n_labels), head(top_down_genes, top_n_labels))
  diffexp$label <- ifelse(diffexp$gene %in% label_genes_df$gene, as.character(diffexp$gene), "")
  
  # Construct legend labels.
  legend_labels <- c("Un-sigif" = "Not Significant",
                     "Down" = paste0("Enriched in ", neg_group_name),
                     "Up" = paste0("Enriched in ", pos_group_name))
  
  # Draw the plot.
  p <- ggplot(diffexp, aes(x = avg_log2FC, y = logP, color = change)) +
    # Adjust point size and opacity by significance status.
    geom_point(aes(size = change, alpha = change)) +
    scale_size_manual(values = c("Un-sigif" = 0.5, "Down" = 1.5, "Up" = 1.5)) +
    scale_alpha_manual(values = c("Un-sigif" = 0.5, "Down" = 0.9, "Up" = 0.9)) +
    scale_color_manual(values = colors, labels = legend_labels) +
    
    # Add dashed threshold lines.
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.7) +
    geom_hline(yintercept = -log10(p_val_adj_threshold), linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.7) +
    
    # Axis labels.
    labs(x = expression(Log[2](Fold~Change)),
         y = expression(-Log[10](Adjusted~P~value)),
         title = title,
         subtitle = paste0("Thresholds: |log2FC| >= ", log2fc_threshold, " & adj.P < ", p_val_adj_threshold)) +
    
    # Style gene labels in italics and avoid overlaps using connector segments.
    geom_text_repel(aes(label = label),
                    size = 4,
                    fontface = "italic",
                    box.padding = unit(0.6, "lines"),
                    point.padding = unit(0.4, "lines"),
                    max.overlaps = 50,
                    segment.color = 'grey50',
                    segment.size = 0.4,
                    min.segment.length = 0,
                    show.legend = FALSE) +
    
    # Add top and bottom padding.
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    
    # Configure the plot theme.
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      axis.line = element_blank() # Use panel.border instead of separate axis lines.
    )
  
  return(p)
}
