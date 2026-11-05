#!/usr/bin/env Rscript
# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Final CSF figures from saved objects or source tables.
# Original source: 10b_revision_myeloid_module_score_compare_layout.R
# Data cohort: csf. Inputs use the compact English data layout.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.


# Replot Myeloid functional module score comparison from existing source tables.

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
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

format_p <- function(p) {
  ifelse(is.na(p), NA_character_, ifelse(p < 0.001, "<0.001", signif(p, 2)))
}

theme_pub <- function(base_size = 9.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black", face = "bold"),
      strip.background = element_rect(fill = "grey94", colour = NA),
      strip.text = element_text(colour = "black", face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.key.size = unit(0.38, "cm"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

save_plot_safe <- function(plot, filename, width, height) {
  pdf_path <- file.path(figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(figure_dir, paste0(filename, ".png"))
  tmp_pdf <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".pdf"))
  tmp_png <- file.path(tempdir(), paste0(filename, "_", Sys.getpid(), ".png"))
  ggsave(tmp_pdf, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(tmp_png, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
  file.copy(tmp_png, png_path, overwrite = TRUE)
  unlink(c(tmp_pdf, tmp_png))
}

# Read the saved object or source table used by this analysis.
sample_qc <- readRDS(file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "CSF_scRNA.rds"))@meta.data %>%
  rownames_to_column("cell") %>%
  count(HLA_Splitter, Group, name = "n_cells") %>%
  rename(Sample = HLA_Splitter) %>%
  arrange(desc(n_cells))
# Export the analysis source data for figure reproduction.
write.csv(sample_qc, file.path(result_dir, "CSF_sample_HLA_Splitter_cell_count_qc.csv"), row.names = FALSE)
valid_samples <- sample_qc$Sample[sample_qc$n_cells > 100]

module_score_sample <- read.csv(
  file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "Myeloid_revision_functional_module_sample_scores.csv"),
  check.names = FALSE
) %>%
  filter(Sample %in% valid_samples) %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    module = factor(module, levels = unique(module))
  )

module_labels <- c(
  "Microglia/phagocytosis" = "Microglia/\nphagocytosis",
  "Complement" = "Complement",
  "Inflammation/chemokine" = "Inflammation/\nchemokine",
  "Antigen presentation" = "Antigen\npresentation",
  "Monocyte activation" = "Monocyte\nactivation",
  "Repair/M2-like" = "Repair/M2-like"
)

module_wilcox_plot <- module_score_sample %>%
  group_by(module) %>%
  group_modify(function(.x, .y) {
    comparisons <- t(utils::combn(group_levels, 2))
    data.frame(
      group1 = comparisons[, 1],
      group2 = comparisons[, 2],
      x = match(comparisons[, 1], group_levels),
      xend = match(comparisons[, 2], group_levels),
      p.value = apply(comparisons, 1, function(pair) {
        values_1 <- .x$module_score[.x$Group == pair[1]]
        values_2 <- .x$module_score[.x$Group == pair[2]]
        if (length(values_1) < 2 || length(values_2) < 2 || length(unique(c(values_1, values_2))) < 2) return(NA_real_)
        suppressWarnings(stats::wilcox.test(values_1, values_2, exact = FALSE)$p.value)
      })
    )
  }) %>%
  ungroup() %>%
  filter(!is.na(p.value), p.value < 0.1) %>%
  left_join(
    module_score_sample %>%
      group_by(module) %>%
      summarise(y_max = max(module_score, na.rm = TRUE), y_min = min(module_score, na.rm = TRUE), .groups = "drop"),
    by = "module"
  ) %>%
  group_by(module) %>%
  arrange(p.value, .by_group = TRUE) %>%
  mutate(
    module = factor(module, levels = levels(module_score_sample$module)),
    y_range = pmax(y_max - y_min, 0.08),
    y.position = y_max + row_number() * y_range * 0.18,
    tick.height = y_range * 0.035,
    label.y = y.position + y_range * 0.065,
    label = format_p(p.value)
  ) %>%
  ungroup()

write.csv(
  module_wilcox_plot,
  file.path(result_dir, "Myeloid_revision_functional_module_wilcox_tests.csv"),
  row.names = FALSE
)

p_module_score <- ggplot(module_score_sample, aes(Group, module_score, fill = Group, colour = Group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, colour = "black", linewidth = 0.42, alpha = 0.34) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 11),
    shape = 21,
    stroke = 0.52,
    size = 2.10,
    alpha = 0.92
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = x, xend = xend, y = y.position, yend = y.position),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = x, xend = x, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_segment(
    data = module_wilcox_plot,
    aes(x = xend, xend = xend, y = y.position, yend = y.position - tick.height),
    inherit.aes = FALSE,
    linewidth = 0.35,
    colour = "black"
  ) +
  geom_text(
    data = module_wilcox_plot,
    aes(x = (x + xend) / 2, y = label.y, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    size = 3.0,
    vjust = -0.30,
    colour = "black"
  ) +
  facet_wrap(~module, scales = "free_y", ncol = 6, labeller = labeller(module = module_labels)) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Module score per sample") +
  theme_pub(10.2) +
  theme(
    axis.title.y = element_text(size = 11.8, face = "bold"),
    axis.text.x = element_text(size = 9.2, colour = "black", angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9.8, colour = "black"),
    strip.text = element_text(size = 9.3, face = "bold", lineheight = 0.95),
    panel.spacing.x = unit(0.75, "lines"),
    legend.position = "none"
  )

save_plot_safe(p_module_score, "Myeloid_functional_module_score_compare", width = 11.0, height = 3.4)

message("Updated Myeloid_functional_module_score_compare with 6 columns x 1 row layout.")
