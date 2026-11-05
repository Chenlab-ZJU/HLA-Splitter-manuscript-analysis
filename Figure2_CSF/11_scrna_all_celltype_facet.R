#!/usr/bin/env Rscript
# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 07_revision_scRNA_all_celltype_facet.R
# Data cohort: csf. Inputs use the compact English data layout.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.


# Restore the original CSF all-celltype UMAP facet plot.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

project_dir <- "D:/项目/2026.Nature Communications/Revison/Data/csf"
figure_dir <- file.path(project_dir, "Figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Read the saved object or source table used by this analysis.
scRNA <- readRDS(file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "CSF_scRNA.rds"))

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_levels <- intersect(group_levels, unique(as.character(scRNA$Group)))

celltype_col <- "Major_celltype"
if (!celltype_col %in% colnames(scRNA@meta.data)) {
  stop("Major_celltype is not available in CSF_scRNA.rds.")
}

celltype_levels <- c(
  "CD4+ T",
  "CD8+ T",
  "Gamma-delta T",
  "B cell",
  "Monocyte",
  "Macrophage",
  "cDC",
  "pDC",
  "Neutrophil",
  "Cycling cell"
)
celltype_levels <- intersect(celltype_levels, unique(as.character(scRNA@meta.data[[celltype_col]])))

celltype_colors <- c(
  "CD4+ T" = "#E42320",
  "CD8+ T" = "#699ECA",
  "Gamma-delta T" = "#731A73",
  "B cell" = "#FFCB5B",
  "Monocyte" = "#0076B9",
  "Macrophage" = "#E89B1E",
  "cDC" = "#0098B2",
  "pDC" = "#8A7355",
  "Neutrophil" = "#8FB4DC",
  "Cycling cell" = "#D65190"
)
celltype_colors <- celltype_colors[celltype_levels]
scRNA$Group <- factor(as.character(scRNA$Group), levels = group_levels)
scRNA$plot_celltype <- factor(as.character(scRNA@meta.data[[celltype_col]]), levels = celltype_levels)

umap_df <- as.data.frame(Embeddings(scRNA, "umap")) %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  mutate(
    Group = scRNA$Group,
    Celltype = scRNA$plot_celltype
  ) %>%
  filter(!is.na(Group), !is.na(Celltype))

x_range <- range(umap_df$UMAP_1, finite = TRUE)
y_range <- range(umap_df$UMAP_2, finite = TRUE)

theme_original_umap <- function() {
  theme_void(base_size = 9, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11, colour = "black"),
      legend.title = element_text(face = "bold", size = 10, colour = "black"),
      legend.text = element_text(size = 8.8, colour = "black"),
      legend.key.size = unit(0.38, "cm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

plot_one_group <- function(group_name) {
  current <- umap_df %>% filter(Group == group_name)
  group_title <- paste0(group_name, " (n=", format(nrow(current), big.mark = ","), ")")
  ggplot() +
    geom_point(
      data = umap_df,
      aes(UMAP_1, UMAP_2),
      colour = "grey72",
      size = 0.14,
      alpha = 0.52,
      stroke = 0
    ) +
    geom_point(
      data = current,
      aes(UMAP_1, UMAP_2, colour = Celltype),
      size = 0.24,
      alpha = 0.82,
      stroke = 0
    ) +
    scale_colour_manual(values = celltype_colors, drop = FALSE, name = "Major cell type") +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    coord_fixed(xlim = x_range, ylim = y_range, expand = FALSE) +
    labs(title = group_title) +
    theme_original_umap()
}

combined_plot <- wrap_plots(lapply(group_levels, plot_one_group), nrow = 1, guides = "collect") &
  theme(legend.position = "right")

save_plot_safe <- function(plot, filename, width, height) {
  target_pdf <- file.path(figure_dir, paste0(filename, ".pdf"))
  target_png <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  file.copy(tmp_pdf, target_pdf, overwrite = TRUE)
  file.copy(tmp_png, target_png, overwrite = TRUE)
  unlink(c(tmp_pdf, tmp_png))
}

save_plot_safe(combined_plot, "scRNA_all_celltype_facet", width = 8.0, height = 5.8)

message("Restored original CSF all-celltype UMAP facet plot: ", file.path(figure_dir, "scRNA_all_celltype_facet.pdf"))
