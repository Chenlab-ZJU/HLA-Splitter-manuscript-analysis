#!/usr/bin/env Rscript
# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 09b_revision_myeloid_umap_cluster_merge_bc_legend.R
# Data cohort: csf. Inputs use the compact English data layout.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.


# Replot Myeloid_UMAP_Cluster with b/c panels merged and legends collected.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(gtools)
  library(patchwork)
  library(ggrepel)
})

project_dir <- "D:/项目/2026.Nature Communications/Revison/Data/csf"

figure_dir <- file.path(project_dir, "Figures")
result_dir <- file.path(project_dir, "Result")

group_levels <- c("Metastases", "Hemorrhage", "Encephalitis")
group_colors <- c(
  "Metastases" = "#E42320",
  "Hemorrhage" = "#699ECA",
  "Encephalitis" = "#4DAF4A"
)
myeloid_colors <- c(
  "Mono" = "#0076B9",
  "HMOX1-Macro" = "#E89B1E",
  "MT1G-Macro" = "#699ECA",
  "PPARG-Macro" = "#4DAF4A",
  "Inflam-Macro" = "#D65190",
  "microglia-like" = "#8A7355",
  "IFN-Myeloid" = "#731A73",
  "Other" = "grey55"
)

theme_pub <- function(base_size = 9.8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black", face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, colour = "black"),
      legend.title = element_text(face = "bold", colour = "black"),
      legend.text = element_text(colour = "black"),
      legend.key.size = unit(0.38, "cm"),
      plot.margin = margin(3, 3, 3, 3)
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

message("Loading existing myeloid object...")
# Read the saved object or source table used by this analysis.
Myeloid <- readRDS(file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "CSF_Myeloid_revision.rds"))

umap_mat <- Embeddings(Myeloid, "umap")
umap_df <- as.data.frame(umap_mat) %>%
  setNames(c("UMAP_1", "UMAP_2")) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    Cluster = factor(Myeloid$seurat_clusters, levels = mixedsort(unique(as.character(Myeloid$seurat_clusters)))),
    Myeloid_celltype = factor(as.character(Myeloid$Myeloid_celltype_revision), levels = names(myeloid_colors)),
    Group = factor(as.character(Myeloid$Group), levels = group_levels)
  ) %>%
  filter(!is.na(Cluster), !is.na(Myeloid_celltype), !is.na(Group))

cluster_centers <- umap_df %>%
  group_by(Cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

cluster_palette <- setNames(
  rep(unname(myeloid_colors), length.out = length(levels(umap_df$Cluster))),
  levels(umap_df$Cluster)
)

p_cluster <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Cluster)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  geom_label_repel(
    data = cluster_centers,
    aes(label = Cluster),
    size = 3.0,
    family = "Arial",
    fontface = "bold",
    label.size = 0.15,
    label.padding = unit(0.10, "lines"),
    fill = scales::alpha("white", 0.82),
    colour = "black",
    seed = 20260813,
    max.overlaps = Inf
  ) +
  scale_colour_manual(values = cluster_palette) +
  labs(title = "Myeloid clusters", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(legend.position = "none", panel.grid = element_blank())

p_anno <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Myeloid_celltype)) +
  geom_point(size = 0.18, alpha = 0.82, stroke = 0) +
  scale_colour_manual(values = myeloid_colors, drop = TRUE, name = "Myeloid state") +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1), ncol = 1)) +
  labs(title = "Myeloid states", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(panel.grid = element_blank())

p_group <- ggplot(umap_df, aes(UMAP_1, UMAP_2, colour = Group)) +
  geom_point(size = 0.18, alpha = 0.78, stroke = 0) +
  scale_colour_manual(values = group_colors, drop = FALSE, name = "Clinical group") +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1), ncol = 1)) +
  labs(title = "Clinical groups", x = "UMAP 1", y = "UMAP 2") +
  theme_pub(9.8) +
  theme(panel.grid = element_blank())

bc_panels <- p_anno + p_group + plot_layout(ncol = 2, guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = unit(0.08, "cm")
  )

combined <- p_cluster + bc_panels + plot_layout(widths = c(1, 2.25))

save_plot_safe(combined, "Myeloid_UMAP_Cluster", width = 10.4, height = 3.65)

message("Updated Myeloid_UMAP_Cluster with merged b/c panels and a single legend column.")
