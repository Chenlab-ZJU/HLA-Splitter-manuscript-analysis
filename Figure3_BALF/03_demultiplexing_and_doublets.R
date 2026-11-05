# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 02_demultiplexing_and_doublets.R
# Data cohort: balf. Inputs use the compact English data layout.

# ==============================================================================
# File: 02_demultiplexing_and_doublets.R
# Description: integrate HLA-Splitter assignments and remove doublets.
# ==============================================================================
rm(list = ls())
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure3_BALF/01_setup_and_helpers.R"))

# Read the saved object or source table used by this analysis.
scRNA_Intergrate <- readRDS("./Result/lavage_scRNA_Intergrate.rds")
filenames <- c("Batch1","Batch2","Batch3","Batch4")

# Integrate HLA-Splitter assignments and draw t-SNE.
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure1_benchmark/01_hla_demultiplex.R"))
sampleInfo <- read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "balf", "HLA-genotypes.csv"), header = TRUE, row.names = 1)
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
SampleInfo <- read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "balf", "SampleInfo.csv"))
scRNA_Intergrate$Group <- SampleInfo$Group[match(scRNA_Intergrate$HLA_Splitter, SampleInfo$Patient.ID)]
scRNA_Intergrate$Group[is.na(scRNA_Intergrate$Group)] <- 'Undefined'
Group_levels <- c("CAP", "LC", "ILD", "Others", "Undefiend")
scRNA_Intergrate$Group <- factor(scRNA_Intergrate$Group, levels = Group_levels)
scRNA_Intergrate$HLA_Splitter = factor(scRNA_Intergrate$HLA_Splitter, 
                                       levels = mixedsort(unique(scRNA_Intergrate$HLA_Splitter)))

# Export the analysis source data for figure reproduction.
saveRDS(scRNA_Intergrate, "Result/lavage_scRNA_Intergrate.rds")

# Draw the overall HLA t-SNE plot.
scRNA_Intergrate = readRDS("Result/lavage_scRNA_Intergrate.rds")
scRNA_Intergrate = scRNA_Intergrate[, scRNA_Intergrate$HLA_Splitter != "Undefined"]

results <- plot_HLA_tSNE(seurat_obj = scRNA_Intergrate, 
                         hla_dir = "./Data/HLA_counts/",
                         filenames = c("Batch1","Batch2","Batch3","Batch4"),
                         out_dir = "./Result/Demultiplex")

df <- data.frame(HLA_demutliplex = scRNA_Intergrate$HLA_Splitter)
df <- table(df) %>% prop.table(1)*100
write.csv(df, paste0("./Result/Demultiplex/demultiplex_res.csv"))


# Exclude cells inconsistent with the configured sex-associated expression criteria.
cell_sample_ids <- scRNA_Intergrate$HLA_Splitter
scRNA_Intergrate$Gender <- Sample_Info$Gender[match(cell_sample_ids, Sample_Info$Patient.ID)]

female_gene <- "XIST"
male_genes <- c("RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "UTY")
target_genes <- intersect(c(female_gene, male_genes), rownames(scRNA_Intergrate))

expr_df <- FetchData(scRNA_Intergrate, vars = target_genes, layer = "data")

# Calculate XIST expression and the mean Y-linked gene score per cell.
if ("XIST" %in% target_genes) {
  xist_expr <- expr_df[["XIST"]]
} else {
  xist_expr <- rep(0, nrow(expr_df))
}

male_genes_present <- intersect(male_genes, target_genes)
if (length(male_genes_present) > 0) {
  y_expr <- rowMeans(expr_df[, male_genes_present, drop = FALSE])
} else {
  y_expr <- rep(0, nrow(expr_df))
}

# Apply the configured sex-associated expression filter.
threshold <- 0.2

is_bad_cell <- rep(FALSE, ncol(scRNA_Intergrate))
for (i in 1:nrow(expr_df)) {
  annotated_gender <- scRNA_Intergrate$Gender[i]
  if (is.na(annotated_gender)) next 
  if (annotated_gender %in% c("Female", "F")) {
    if (y_expr[i] > threshold) is_bad_cell[i] <- TRUE
  } else if (annotated_gender %in% c("Male", "M")) {
    if (xist_expr[i] > threshold) is_bad_cell[i] <- TRUE
  }
}

bad_cells_names <- rownames(expr_df)[is_bad_cell]
message(paste("Removed", length(bad_cells_names), "cells due to mismatched gender gene expression."))

scRNA_Cleaned <- subset(scRNA_Intergrate, cells = setdiff(colnames(scRNA_Intergrate), bad_cells_names))
scRNA_Cleaned <- NormalizeData(scRNA_Cleaned) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()

# Compare donor assignments with sex-associated gene expression.
gene_matrix <- AverageExpression(scRNA_Cleaned, group.by = c("HLA_Splitter"), assays = "RNA", return.seurat = TRUE)
gene_matrix <- gene_matrix@assays$RNA$scale.data %>% as.matrix()

Sample_Info = read.csv(file.path("D:/项目/2026.Nature Communications/Revison/Data", "balf", "SampleInfo.csv"))
annotation_col = data.frame(Gender = Sample_Info$Gender[match(colnames(gene_matrix), Sample_Info$Patient.ID)])

cell_counts <- table(scRNA_Cleaned$HLA_Splitter)
annotation_col$Cell_num = as.numeric(cell_counts[colnames(gene_matrix)])
rownames(annotation_col) = colnames(gene_matrix)

# XIST is associated with X inactivation; Y-linked genes provide complementary sex-associated expression checks.
gender_genes <- c("XIST", "RPS4Y1", "EIF1AY", "DDX3Y", "KDM5D", "UTY")

gender_genes <- intersect(gender_genes, rownames(gene_matrix))

sample_order <- order(annotation_col$Gender)
mat <- gene_matrix[gender_genes, sample_order, drop = FALSE]
annotation_col <- annotation_col[sample_order, , drop = FALSE]

plot_heatmap = pheatmap(mat, 
                        cluster_rows = TRUE,    # Cluster rows to organize X- and Y-linked gene patterns.
                        cluster_cols = FALSE,   # Disable column clustering to retain the specified donor order.
                        border_color = "black", 
                        fontsize = 10, 
                        cellwidth = 20, 
                        cellheight = 20,
                        annotation_col = annotation_col)

ggsave(paste0("./Result/demultiplex/Gender_related_GENEs.pdf"), plot_heatmap, dpi = 300, height = 6, width = 15)

