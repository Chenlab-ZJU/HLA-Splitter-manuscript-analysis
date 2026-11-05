# Purpose: Transitive shared function dependency.
# Original source: plot_dotplot.R
# Data cohort: shared. Inputs use the compact English data layout.

#' Create a dot plot for features expression in scRNA-seq data.
#'
#' This function featuresrates a dot plot to visualize the average expression and
#' percentage of cells expressing a features of interest across different groups
#' and cell subtypes from a Seurat object.
#'
#' @param seurat_object A Seurat object containing scRNA-seq data.
#' @param features A character string specifying the features to plot.
#' @param group_by A character string specifying the metadata column for
#'   grouping the cells (e.g., "Group"). Defaults to "Group".
#' @param cluster_by A character string specifying the metadata column for
#'   cell subtypes (e.g., "revised_celltype"). Defaults to "revised_celltype".
#' @param color_low A character string for the low end of the color gradient.
#'   Defaults to "lightblue".
#' @param color_high A character string for the high end of the color gradient.
#'   Defaults to "red".
#'
#' @return A ggplot2 object representing the dot plot.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom dplyr %>% mutate group_by across all_of summarise select
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient
#' @importFrom ggplot2 scale_size_continuous labs theme_classic theme
#' @importFrom ggplot2 element_text
#' @importFrom tidyr drop_na
#' @export
#'
plot_dotplot <- function(seurat_object,
                         features,
                         group_by = NULL,
                         cluster_by = NULL
                         ) {
  
  # Prepare data and calculate summaries.
  # Check if the feature(s) exist as genes in the Seurat object.
  is_gene_feature <- base::all(features %in% base::rownames(seurat_object))
  # Check if the feature(s) exist in the object's metadata.
  is_meta_feature <- base::all(features %in% base::colnames(seurat_object@meta.data))
  
  if (is_gene_feature) {
    # If the feature is a gene, extract the normalized expression data.
    expression_data <- Seurat::GetAssayData(seurat_object, assay = "RNA", layer = "data")[features, , drop = FALSE]
    base::cat("Extracting gene expression data.\n")
    
  } else if (is_meta_feature) {
    if (!base::is.numeric(seurat_object@meta.data[[features]])) {
      base::stop("The metadata feature is not a numeric value. It cannot be processed for quantitative analysis.")
    }
    expression_data <- base::as.data.frame(seurat_object@meta.data[, features, drop = FALSE])
    base::cat("Extracting numerical metadata feature(s).\n")
    
  } else {
    # Stop the function if the feature is not found in either genes or metadata.
    base::stop("Feature(s) not found in either the Seurat object's gene expression data or its metadata.")
  }
  
  # Extract metadata.
  meta_data <- seurat_object@meta.data
  expression_vector <- base::as.numeric(unlist(expression_data))
  
  # Join expression values with metadata.
  df <- meta_data %>%
    dplyr::mutate(featuresExpression = expression_vector)
  
  # Calculate mean expression and positive-cell fraction for each subtype and group.
  dotplot_data <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_by, cluster_by)))) %>%
    dplyr::summarise(
      AvgExpression = base::mean(featuresExpression),
      PctExpressed = base::sum(featuresExpression > 0) / dplyr::n() * 100,
      .groups = "drop"
    ) %>%
    # Remove rows with missing AvgExpression.
    tidyr::drop_na(AvgExpression)
  
  # Draw the dot plot.
  p = ggplot2::ggplot(dotplot_data, ggplot2::aes(x = .data[[cluster_by]], y = .data[[group_by]])) +
    ggplot2::geom_point(ggplot2::aes(size = PctExpressed, color = AvgExpression)) +
    
    # Configure the color scale.
    ggplot2::scale_color_gradientn(colours = c("navy","white","firebrick3")) +
    
    # Set the point-size range.
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    
    # Add titles and labels.
    ggplot2::labs(
      title = base::paste0(features, " Expression Dot Plot"),
      x = "Cell Subtype",
      y = "Group",
      size = "Cell percent (%)",
      color = "Average Expression"
    ) +
    
    # Adjust the plot theme and layout.
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold")
    )
  
  return(p)
}
