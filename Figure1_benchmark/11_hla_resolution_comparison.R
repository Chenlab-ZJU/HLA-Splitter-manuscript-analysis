# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream HLA typing, resolution and PBMC demultiplexing.
# Original source: Compare_HLA_resolution.R
# Data cohort: pbmc30. Inputs use the compact English data layout.

setwd(file.path("D:/项目/2026.Nature Communications/Revison/Data", "pbmc30"))
rm(list = ls())

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(reshape2)
library(R.utils)
library(pheatmap)
library(ggridges)
library(stringr)
library(ggsci)
library(gtools)
library(scales)
library(patchwork)

# Define the ggsci IGV palette for donor labels.
color <- pal_igv("default")(51)

filename = "30PBMC"

my_theme <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = "black"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Read HLA-based donor assignments.
# Read the saved object or source table used by this analysis.
demultiplex_res_1_field = read.csv("./Result/HLA_Splitter_1field/HLA_demultiplex.csv")
demultiplex_res_2_field = read.csv("./Result/HLA_Splitter_2field/HLA_demultiplex.csv")

# ==============================================================================
# Compare assignments using 1-field and 2-field HLA references.
# ==============================================================================
comparison_df <- inner_join(
  demultiplex_res_1_field %>% select(Cell.barcode, prediction_1field = Sample_prediction),
  demultiplex_res_2_field %>% select(Cell.barcode, prediction_2field = Sample_prediction),
  by = "Cell.barcode"
)

# Calculate assignment concordance.
comparison_df <- comparison_df %>%
  mutate(is_consistent = (prediction_1field == prediction_2field))

# Summarize overall agreement.
consistency_summary <- comparison_df %>%
  summarise(
    Total_Cells = n(),
    Consistent_Cells = sum(is_consistent, na.rm = TRUE),
    Consistency_Rate = Consistent_Cells / Total_Cells
  )

print("=== Assignment agreement between 1-field and 2-field references ===")
print(consistency_summary)

# ==============================================================================
# Plot assignments from the two HLA resolutions as bars.
# ==============================================================================
counts_1f <- comparison_df %>%
  dplyr::count(prediction_1field, name = "Count") %>%
  mutate(Method = "1-field", Sample = prediction_1field) %>%
  select(Sample, Method, Count)

counts_2f <- comparison_df %>%
  dplyr::count(prediction_2field, name = "Count") %>%
  mutate(Method = "2-field", Sample = prediction_2field) %>%
  select(Sample, Method, Count)

df_counts <- bind_rows(counts_1f, counts_2f)

df_consistency <- comparison_df %>%
  filter(!grepl("Doublet|Negative", prediction_1field, ignore.case = TRUE)) %>%
  group_by(Sample = prediction_1field) %>%
  summarise(
    Total = n(),
    Consistent = sum(is_consistent, na.rm = TRUE),
    Consistency_Rate = Consistent / Total
  )

# Use the same sample factor order on the x axis.
sample_levels <- mixedsort(unique(df_counts$Sample))
df_counts$Sample <- factor(df_counts$Sample, levels = sample_levels)
df_consistency$Sample <- factor(df_consistency$Sample, levels = sample_levels)

p_counts <- ggplot(df_counts, aes(x = Sample, y = Count, fill = Method)) +
  geom_col(position = position_dodge(width = 0.8), color = "black", alpha = 0.85) +
  scale_fill_manual(values = c("1-field" = "#4DAF4A", "2-field" = "#699ECA")) +
  labs(
    title = "Cell Counts per Sample (1-field vs 2-field)",
    y = "Number of Cells"
  )+ my_theme()

p_rate <- ggplot(df_consistency, aes(x = Sample, y = Consistency_Rate)) +
  geom_col(fill = "#E42320", color = "black", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = percent(Consistency_Rate, accuracy = 0.1)), 
            vjust = -0.5, size = 3.5) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1.1)) + # Leave y-axis space for labels.
  labs(
    title = "Demultiplexing Consistency Rate per Sample",
    x = "Predicted Sample (Baseline: 1-field)",
    y = "Consistency Rate"
  ) + my_theme()

combined_plot <- p_counts / p_rate + 
  plot_layout(heights = c(1, 1))

print(combined_plot)

ggsave("./Figures/Demultiplex_Counts_and_Consistency.pdf", plot = combined_plot, width = 14, height = 6)



