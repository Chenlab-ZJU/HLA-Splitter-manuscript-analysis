#!/usr/bin/env Rscript
# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 10c_revision_myeloid_heatmap_label.R
# Data cohort: csf. Inputs use the compact English data layout.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.


# Replot Myeloid group DEG key-gene heatmap from existing source table.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

project_dir <- "D:/项目/2026.Nature Communications/Revison/Data/csf"

figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")

theme_pub <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black", face = "bold"),
      plot.title = element_text(colour = "black", face = "bold"),
      strip.text = element_text(colour = "black", face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

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

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")

# Read the saved object or source table used by this analysis.
heat_df <- read.csv(
  file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "Myeloid_revision_group_DEG_heatmap_source.csv"),
  check.names = FALSE
) %>%
  mutate(
    module = factor(module, levels = unique(module)),
    module_gene = factor(module_gene, levels = rev(unique(module_gene))),
    Group = factor(Group, levels = group_levels),
    z_score = as.numeric(z_score)
  )

p_heat <- ggplot(heat_df, aes(Group, module_gene, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#6F99C3", mid = "#F7F7F7", high = "#B7192B", midpoint = 0, name = "Z-score") +
  scale_y_discrete(labels = function(x) sub("^.*__", "", x)) +
  labs(title = "Myeloid cells", x = NULL, y = NULL) +
  theme_pub(9.2) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 10.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 7.6),
    strip.text.y = element_text(size = 7.6, angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey93", colour = NA),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot_safe(p_heat, "Myeloid_group_DEG_key_gene_heatmap", width = 4.7, height = 7.05)

message("Updated Myeloid_group_DEG_key_gene_heatmap with Myeloid cells label.")
