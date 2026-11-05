# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: 01_data_preprocessing_demux.R
# Data cohort: balf. Inputs use the compact English data layout.

# ==============================================================================
# File: 01_data_preprocessing_demux.R
# Description: merge data, perform QC, visualize assignments and remove doublets with Scrublet.
# ==============================================================================

rm(list = ls())
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure3_BALF/01_setup_and_helpers.R"))
dir.create("Figures/QC", recursive = TRUE, showWarnings = FALSE)

# Read the single-cell expression matrix and review QC.
filenames <- c("Group1","Group2","Group3","Group4")
scRNA_list <- lapply(seq_along(filenames), function(i) {
  Counts <- Read10X(data.dir = paste0("./Data/GEM/", filenames[i], "/filtered_feature_bc_matrix/"))
  obj <- CreateSeuratObject(Counts, project = filenames[i], assay = "RNA")
  obj <- RenameCells(obj, new.names = paste0(str_sub(colnames(obj), 1, 16), "-", i))
  return(obj)
})
scRNA_Intergrate <- merge(scRNA_list[[1]], y = scRNA_list[-1], project = "lavage_Total")


# Perform single-cell QC and filtering.
scRNA_Intergrate[["percent.mt"]] <- PercentageFeatureSet(scRNA_Intergrate, pattern = "^MT-")
scRNA_Intergrate <- subset(scRNA_Intergrate, subset = nCount_RNA > 200 & percent.mt < 20)

# Plot quality-control violins.
qc_metrics <- list(nFeature = c("nFeature_RNA", 5000), 
                   nCount = c("nCount_RNA", 10000), 
                   MT = c("percent.mt", 20))

dir.create("./Figures/QC", recursive = TRUE, showWarnings = FALSE)
for(name in names(qc_metrics)) {
  p <- ggplot(scRNA_Intergrate@meta.data, aes(x = orig.ident, y = .data[[qc_metrics[[name]][1]]])) +
    geom_violin(aes(fill = orig.ident), width = 1, alpha = 0.5) +
    geom_jitter(aes(color = orig.ident), width = 0.2, size = 0.01) +
    geom_boxplot(fill = "white", width = 0.3, outlier.size = 0) +
    ylim(c(0, as.numeric(qc_metrics[[name]][2]))) + theme_bw()
  ggsave(paste0("./Figures/QC/Vln_", name, ".pdf"), plot = p, height = 5, width = 8)
}

# Export the analysis source data for figure reproduction.
saveRDS(scRNA_Intergrate, "./Result/lavage_scRNA_Intergrate.rds")
