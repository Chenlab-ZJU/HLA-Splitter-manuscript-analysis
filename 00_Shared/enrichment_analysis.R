# Purpose: Transitive shared function dependency.
# Original source: Enrichment_analysis.R
# Data cohort: shared. Inputs use the compact English data layout.

#' Perform Gene Enrichment Analysis for a Two-Group Comparison
#'
#' This function takes a data.frame of differentially expressed genes from a two-group
#' comparison (e.g., Seurat FindMarkers), performs enrichment analysis for the up-regulated
#' and down-regulated genes separately, and generates a faceted bar plot.
#'
#' @param diffexp A data.frame of differential expression results containing 'p_val_adj', 'avg_log2FC'.
#' @param pos_group_name Name of the group with positive log2FC (Up-regulated).
#' @param neg_group_name Name of the group with negative log2FC (Down-regulated).
#' @param dbs A character vector of database names for enrichR.
#' @param top_n_pathways The number of top enriched pathways to display per group (default: 6).
#' @param plot_colors Colors for the positive and negative groups.
#' @param p_val_adj_threshold The adjusted p-value threshold (default: 0.05).
#' @param log2fc_threshold The log2 fold change threshold (default: 1).
#' @param output_csv_path File path to save the full enrichment results.
#' @return A ggplot object of the enrichment results.
#' @importFrom enrichR enrichr
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual coord_flip labs theme_minimal theme element_blank element_text element_line margin facet_wrap
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_replace_all

Enrichment_analysis<- function(diffexp,
                               pos_group_name = "Group A",
                               neg_group_name = "Group B",
                               dbs = c("GO_Biological_Process_2023"),
                               top_n_pathways = 6, 
                               plot_colors = c("#e74c3c", "#3498db"),
                               p_val_adj_threshold = 0.05,
                               log2fc_threshold = 1,
                               output_csv_path = "./TwoGroup_Enrichment_Results.csv") {
  
  if (!requireNamespace("enrichR", quietly = TRUE)) stop("Package 'enrichR' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  
  # Ensure a gene column exists when FindMarkers stores symbols as row names.
  if (!"gene" %in% colnames(diffexp)) {
    diffexp$gene <- rownames(diffexp)
  }
  
  # Extract upregulated and downregulated gene lists.
  genes_up <- diffexp$gene[diffexp$p_val_adj < p_val_adj_threshold & diffexp$avg_log2FC >= log2fc_threshold]
  genes_down <- diffexp$gene[diffexp$p_val_adj < p_val_adj_threshold & diffexp$avg_log2FC <= -log2fc_threshold]
  
  # Run enrichR and format results using an internal helper.
  run_enrichment <- function(gene_list, group_label) {
    if (length(gene_list) == 0) return(NULL)
    
    res <- enrichr(gene_list, dbs)
    
    # Combine results across databases.
    df_list <- lapply(names(res), function(db) {
      if (nrow(res[[db]]) > 0) {
        temp <- res[[db]]
        temp$Database <- db
        return(temp)
      }
      return(NULL)
    })
    
    enriched_df <- do.call(rbind, df_list[!sapply(df_list, is.null)])
    if (is.null(enriched_df) || nrow(enriched_df) == 0) return(NULL)
    
    # Select required columns and calculate -log10(P_value).
    df_processed <- data.frame(
      Term = enriched_df$Term,
      P_value = enriched_df$P.value,
      Genes = enriched_df$Genes,
      Group = group_label,
      Database = enriched_df$Database
    )
    
    df_processed <- df_processed[order(df_processed$P_value, decreasing = FALSE), ]
    df_processed$minus_log_pvalue <- -log10(df_processed$P_value)
    return(df_processed)
  }
  
  # Run enrichment analysis.
  message("Analyzing enriched pathways for: ", pos_group_name)
  res_up <- run_enrichment(genes_up, pos_group_name)
  
  message("Analyzing enriched pathways for: ", neg_group_name)
  res_down <- run_enrichment(genes_down, neg_group_name)
  
  # Combine and save complete results.
  combined_all_results <- bind_rows(res_up, res_down)
  if (nrow(combined_all_results) == 0) {
    warning("No enriched pathways found for either group.")
    return(NULL)
  }
# Export the analysis source data for figure reproduction.
  write.csv(combined_all_results, output_csv_path, row.names = FALSE)
  message("Full enrichment results saved to: ", output_csv_path)
  
  # Plot the top N terms for each group.
  top_results_list <- list()
  if (!is.null(res_up)) top_results_list[[pos_group_name]] <- head(res_up, top_n_pathways)
  if (!is.null(res_down)) top_results_list[[neg_group_name]] <- head(res_down, top_n_pathways)
  
  combined_top_results <- bind_rows(top_results_list)
  
  # Set factor levels to control plotting order.
  combined_top_results$Group <- factor(combined_top_results$Group, levels = c(pos_group_name, neg_group_name))
  combined_top_results$Term <- factor(combined_top_results$Term, levels = rev(unique(combined_top_results$Term)))
  
  # Assign colors.
  names(plot_colors) <- c(pos_group_name, neg_group_name)
  
  # Draw faceted bar plots.
  p <- ggplot(combined_top_results, aes(x = Term, y = minus_log_pvalue, fill = Group)) +
    geom_col(linewidth = 0.5, alpha = 0.7, color = "black") +
    scale_fill_manual(values = plot_colors) +
    coord_flip() +
    labs(
      x = "Enriched Pathway",
      y = expression(-Log[10](P~value)),
      fill = "Enriched in"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_blank(), # Hide y-axis text because pathway labels are placed inside the bars.
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(), 
      strip.text = element_text(size = 14, face = "bold"), # Facet titles.
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1, "lines"),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    ) +
    # Place pathway names inside bars; omit GO identifiers from labels.
    geom_text(aes(label = sub("\\s*\\(.*", "", Term), y = 0.2),
              hjust = 0, 
              color = "black", 
              size = 4.5, 
              fontface = "bold") + 
    # Use free y scales to give each group its own plotting space.
    facet_wrap(~Group, scales = "free_y", ncol = 1, strip.position = "top")
  
  return(p)
}
