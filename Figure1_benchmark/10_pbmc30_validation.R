# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream HLA typing, resolution and PBMC demultiplexing.
# Original source: 30PBMC.R
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

dir.create("Figures", showWarnings = FALSE)
dir.create("Result", showWarnings = FALSE)
dir.create("Result/demultiplex", showWarnings = FALSE)

# Define the ggsci IGV palette for donor labels.
color <- pal_igv("default")(51)

filename = "30PBMC"

# ==========================================
# Read data and perform initial quality control.
# ==========================================
Counts <- Read10X(data.dir = paste0("./Data/GEM/filtered_feature_bc_matrix/"))
scRNA_all <- CreateSeuratObject(Counts, project = filename, assay = "RNA")

scRNA_all[["percent.mt"]] <- PercentageFeatureSet(scRNA_all, pattern = "^MT-")

p1 = VlnPlot(scRNA_all, features=c("nFeature_RNA","nCount_RNA","percent.mt"), ncol=2, pt.size=0)
p2 = FeatureScatter(scRNA_all, feature1="nCount_RNA", feature2="nFeature_RNA")
ggsave(paste0("./Figures/QC.pdf"), p1+p2, dpi=300, height=10, width = 12)

# Filter low-quality cells.
scRNA_all <- subset(scRNA_all, subset = nFeature_RNA > 200 & nFeature_RNA < 100000 & percent.mt < 20)

# ==========================================
# Cluster cells and identify doublets with Scrublet.
# ==========================================
scRNA_all <- NormalizeData(scRNA_all) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
scRNA_all <- RunTSNE(scRNA_all, dims = 1:30)
scRNA_all <- RunUMAP(scRNA_all, dims = 1:30)
scRNA_all <- FindNeighbors(scRNA_all, dims = 1:30) %>% FindClusters(resolution=0.8)

# install.packages("remotes")
# remotes::install_github("r-lib/devtools")
# devtools::install_github("furlan-lab/scrubletR")
library(scrubletR)
scRNA_all <- scrublet(scRNA_all)
p1 = FeaturePlot(scRNA_all, features = "doublet_scores", reduction = "tsne")
ggsave("./Figures/scrublet_doublet_score.pdf", p1, height=5, width = 6)

scRNA_all$doublets <- scRNA_all$doublet_scores > 0.2   
p2 = DimPlot(scRNA_all, group.by = "doublets", cols=c("goldenrod", "navy"), reduction = "tsne")
ggsave("./Figures/scrublet_doublet.pdf", p2, height=5, width = 6)

# ==========================================
# Retain singlets and recluster.
# ==========================================
scRNA = scRNA_all[, !scRNA_all$doublets]
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData() %>% RunPCA()
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- FindNeighbors(scRNA, dims = 1:30) %>% FindClusters(resolution=0.3)

# ==========================================
# Identify and export subtype differential-expression results.
# ==========================================
scRNA <- SetIdent(scRNA, value = "seurat_clusters")
# Discover marker genes with the original assay and comparison settings.
cluster.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# Export the analysis source data for figure reproduction.
write.csv(cluster.markers, "./Result/scRNA_markers.csv", row.names = FALSE)

df_summary <- data.frame(seurat_clusters = levels(scRNA$seurat_clusters))
# Select the top 20 genes per cluster.
df_top_genes <- cluster.markers %>%
  filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) %>% 
  summarise(top_genes = paste0(gene, collapse = ", ")) %>% 
  ungroup() %>% 
  rename(seurat_clusters = cluster)

final_df <- left_join(df_summary, df_top_genes, by = "seurat_clusters")
final_df$Celltype = c("")
write.csv(final_df, "./Result/scRNA_celltype.csv", row.names = FALSE)

message("Open ./Result/scRNA_markers_summary.csv, add the Celltype column, and save before continuing.")

# ==========================================
# Review and rename subtypes using highly expressed markers.
# ==========================================
# Read the saved object or source table used by this analysis.
df = read.csv("./Result/scRNA_celltype.csv")

# Construct the named annotation mapping.
Celltype = df$Celltype
names(Celltype) = as.character(df$seurat_clusters)

# Replace numeric cluster labels with cell-type annotations.
scRNA <- SetIdent(scRNA, value = "seurat_clusters")
scRNA <- RenameIdents(scRNA, Celltype)
scRNA$Celltype <- Idents(scRNA)

# Set factor levels to control legend order.
# scRNA$Celltype <- factor(scRNA$Celltype, levels = c("CD4+ T", "CD8+ T", "B cell", "NK",
#                                           "cDC", "pDC", "Monocyte", "Macrophage",
#                                           "Neutrophil", "Cycling cell", "Low-quality"))

# ==========================================
# Assemble UMAP panels.
# ==========================================
my_colors <- pal_igv("default")(length(levels(scRNA$Celltype)))

# Panel 1: numeric clusters.
scRNA <- SetIdent(scRNA, value = "seurat_clusters")
p1 <- DimPlot(scRNA, reduction = "umap", label = TRUE, alpha = 0.5, label.size = 3) +
  theme(legend.position = "right", legend.box = "vertical", legend.key.size = unit(3, "pt")) +
  scale_color_discrete() 

# Panel 2: annotated cell types.
scRNA <- SetIdent(scRNA, value = "Celltype")
p2 <- DimPlot(scRNA, reduction = "umap", label = TRUE, alpha = 0.5, label.size = 3) +
  theme(legend.position = "right", legend.box = "vertical", legend.key.size = unit(3, "pt")) +
  scale_color_manual(values = my_colors) # Use the predefined palette.

# Assemble and save the figure.
combined_plot <- p1 + p2
ggsave("./Figures/scRNA_UMAP_celltype.pdf", combined_plot, height = 5, width = 12)


# ==========================================
# Demultiplex samples using HLA-Splitter.
# ==========================================
# Read HLA allele-count inputs.
sampleInfo = read.csv("./Data/HLA_list/HLA-genotypes_ONT.csv", header=TRUE, row.names = 1)
namelist = colnames(sampleInfo)

# Construct the donor HLA list.
HLA_list = list()
for (i in 1:length(namelist)){ HLA_list[[i]] = unique(sampleInfo[,i]) }
names(HLA_list) = namelist

# Read the count matrix.
HLA_label = read.delim("./Data/HLA_counts/sorted_genotypes.txt", header = F)
barcodes = read.delim(paste0("./Data/HLA_counts/output.barcodes.txt"), header=F)
HLA_matrix = readMM(paste0("./Data/HLA_counts/HLA_type.mtx")) %>% as.matrix() %>% as.data.frame()
colnames(HLA_matrix) = HLA_label$V1
rownames(HLA_matrix) = paste0(barcodes$V1, "-1")

# Load the demultiplexing function and run assignments.
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure1_benchmark/01_hla_demultiplex.R"))
sample_predict <- HLA_demultiplex(HLA_list = HLA_list, 
                                  HLA_matrix = HLA_matrix, 
                                  graph_corr = TRUE, 
                                  outdir = "Result/demultiplex")

# ==========================================
# Map assignments to the Seurat object and plot results.
# ==========================================
scRNA$HLA_splitter = sample_predict[match(names(sample_predict), colnames(scRNA))] # Match cells by barcode.
scRNA$HLA_splitter[is.na(scRNA$HLA_splitter)] = "Undefined"

# Draw UMAPs for clusters and annotated cell types.
scRNA = SetIdent(scRNA, value = "seurat_clusters")
p1 = DimPlot(scRNA, reduction = "umap", label=T) +
  scale_color_discrete(type = c(color)) +
  theme(legend.position = "right", legend.box = "vertical", legend.key.size = unit(3,"pt"))

scRNA = SetIdent(scRNA, value = "HLA_splitter")
p2 = DimPlot(scRNA, reduction = "umap", label=F) +
  scale_color_discrete(type = c(color)) +
  theme(legend.position = "right", legend.box = "vertical", legend.key.size = unit(3,"pt"))

# Assemble and save panels.
ggsave(paste0("./Figures/UMAP_Cluster_Splitter.pdf"), p1 + p2, height = 5, width = 12)

# ==========================================
# Plot cell-type proportions.
# ==========================================
df_ratio = data.frame(Sample = scRNA$HLA_splitter, Celltype = scRNA$Celltype) %>% 
  table() %>% prop.table(margin = 1) %>% data.frame() # Use margin=1 to calculate proportions within each sample.
df_ratio$Freq = df_ratio$Freq * 100

p3 = ggplot(df_ratio, aes(x = Sample, y = Freq, fill = Celltype)) +
  geom_col(position = "fill", alpha=0.7) +
  scale_fill_manual(values = color) +
  theme_classic() +
  labs(y = "Frequency of cell types (%)", x = "") +
  theme(
    text = element_text(size = 15),
    axis.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12), 
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
  )
ggsave(paste0("./Figures/Celltype_ratio.pdf"), p3, height = 5, width = 10)

# ==========================================
# Save output files.
# ==========================================
saveRDS(scRNA_all, "./Code/30PBMC_all.rds")
saveRDS(scRNA, "./Code/30PBMC.rds")
