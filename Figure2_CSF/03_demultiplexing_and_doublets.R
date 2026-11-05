# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream cohort processing and reusable helper definitions.
# Original source: 02_demultiplexing_and_doublets.R
# Data cohort: csf. Inputs use the compact English data layout.
# Manuscript v3: retain samples with >100 confidently assigned cells.
# Sample-level comparisons are descriptive in v3 (3/2/2 samples).
# Legacy inferential annotations require author reconciliation before release.

# ==============================================================================
# File: 02_demultiplexing_and_doublets.R
# Description: integrate HLA-Splitter assignments and remove doublets.
# ==============================================================================
rm(list = ls())
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure2_CSF/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA_Intergrate <- readRDS("./Result/CSF_scRNA_Intergrate.rds")
filenames <- c("Batch1","Batch2","Batch3")

# Integrate HLA-Splitter donor assignments.
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure1_benchmark/01_hla_demultiplex.R"))
sampleInfo <- read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "HLA-genotypes.csv"), header = TRUE, row.names = 1)
HLA_list <- lapply(sampleInfo, unique)

for (i in seq_along(filenames)) {
  base_path <- paste0("./Data/HLA_counts/", filenames[i])
  HLA_label <- read.delim(paste0(base_path, "/sorted_genotypes.txt"), header = FALSE)
  barcodes <- read.delim(paste0(base_path, "/output.barcodes.txt"), header = FALSE)
  
  HLA_matrix <- readMM(paste0(base_path, "/HLA_type.mtx")) %>% as.matrix() %>% as.data.frame()
  colnames(HLA_matrix) <- HLA_label$V1
  rownames(HLA_matrix) <- barcodes$V1
  
  outdir <- paste0("./Result/Demultiplex/", filenames[i])
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  sample_predict <- HLA_demultiplex(HLA_list = HLA_list, HLA_matrix = HLA_matrix, 
                                    graph_corr = TRUE, outdir = outdir)
}

# Add demultiplexing results to cell metadata.
demultiplex_result <- do.call(rbind, lapply(seq_along(filenames), function(i) {
  tmp <- read.csv(paste0("./Result/Demultiplex/", filenames[i], "/HLA_demultiplex.csv"))
  tmp$Cell.barcode <- paste0(tmp$Cell.barcode, "-", i)
  return(tmp)
}))
rownames(demultiplex_result) <- demultiplex_result$Cell.barcode

scRNA_Intergrate$HLA_Splitter <- demultiplex_result[colnames(scRNA_Intergrate), "Sample_prediction"]
scRNA_Intergrate$HLA_Splitter[is.na(scRNA_Intergrate$HLA_Splitter)] <- "Undefined"

# Remove doublets with scrubletR.
options(future.globals.maxSize = 10000 * 1024^2)

scRNA_Intergrate <- NormalizeData(scRNA_Intergrate) %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData() %>% RunPCA() %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  JoinLayers()

scRNA_Intergrate <- scrublet(scRNA_Intergrate)
scRNA_Intergrate <- subset(scRNA_Intergrate, subset = doublet_scores < 0.2)

# Add clinical group information and save the object.
SampleInfo <- read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "SampleInfo.csv"))
scRNA_Intergrate$Group <- SampleInfo$Group[match(scRNA_Intergrate$HLA_Splitter, SampleInfo$Patient.ID)]
scRNA_Intergrate$Group[is.na(scRNA_Intergrate$Group)] <- 'Undefined'
Group_levels <- c("Metastases", "Hemorrhage", "Encephalitis", "Undefined")
scRNA_Intergrate$Group <- factor(scRNA_Intergrate$Group, levels = Group_levels)
scRNA_Intergrate$HLA_Splitter = factor(scRNA_Intergrate$HLA_Splitter, 
                                       levels = mixedsort(unique(scRNA_Intergrate$HLA_Splitter)))

# Export the analysis source data for figure reproduction.
saveRDS(scRNA_Intergrate, "Result/CSF_scRNA_Intergrate.rds")


# Draw the overall HLA t-SNE plot.
#scRNA_Intergrate = readRDS("Result/CSF_scRNA_Intergrate.rds")
scRNA_Intergrate = scRNA_Intergrate[, scRNA_Intergrate$HLA_Splitter != "Undefined"]

results <- plot_HLA_tSNE(seurat_obj = scRNA_Intergrate, 
                         hla_dir = "./Data/HLA_counts/",
                         filenames = c("Batch1","Batch2","Batch3"),
                         out_dir = "./Result/Demultiplex",
                         base_colors = GLOBAL_COLORS)

df <- data.frame(HLA_demutliplex = scRNA_Intergrate$HLA_Splitter)
df <- table(df) %>% prop.table(1)*100
write.csv(df, paste0("./Result/Demultiplex/demultiplex_res.csv"))

# Compare sex-associated expression with donor assignments.
#scRNA_Intergrate = readRDS("Result/CSF_scRNA_Intergrate.rds")
scRNA_Intergrate = scRNA_Intergrate[, scRNA_Intergrate$HLA_Splitter != "Undefined"]

gene_matrix <- AverageExpression(scRNA_Intergrate, group.by = c("HLA_Splitter"), assays = "RNA", return.seurat = TRUE)
gene_matrix <- gene_matrix@assays$RNA$scale.data %>% as.matrix()

Sample_Info = read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "csf", "SampleInfo.csv"))
annotation_col = data.frame(Gender = Sample_Info$Gender[match(colnames(gene_matrix), Sample_Info$Patient.ID)])
 
cell_counts <- table(scRNA_Intergrate$HLA_Splitter)
annotation_col$Cell_num = as.numeric(cell_counts[colnames(gene_matrix)])
rownames(annotation_col) = colnames(gene_matrix)

# XIST is associated with X inactivation; Y-linked genes provide complementary sex-associated expression checks.
gender_genes <- c("XIST", "RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "UTY")

gender_genes <- intersect(gender_genes, rownames(gene_matrix))

sample_order <- order(annotation_col$Gender)
mat <- gene_matrix[gender_genes, sample_order, drop = FALSE]
annotation_col <- annotation_col[sample_order, , drop = FALSE]

library(pheatmap)
plot_heatmap = pheatmap(mat, 
                        cluster_rows = TRUE,    # Cluster rows to organize X- and Y-linked gene patterns.
                        cluster_cols = FALSE,   # Disable column clustering to retain the specified donor order.
                        border_color = "black", 
                        fontsize = 10, 
                        cellwidth = 20, 
                        cellheight = 20,
                        annotation_col = annotation_col)

ggsave(paste0("./Result/demultiplex/Gender_related_GENEs.pdf"), plot_heatmap, dpi = 300, height = 6, width = 10)
