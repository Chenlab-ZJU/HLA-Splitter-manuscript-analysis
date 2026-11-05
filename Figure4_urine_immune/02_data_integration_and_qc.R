# Edit input, output and helper-script paths below for your environment.
# Run the required sections manually; no shared execution entry point is required.
# Purpose: Upstream urine processing and sample-level signature analysis.
# Original source: 01_data_integration_and_QC.R
# Data cohort: urine. Inputs use the compact English data layout.

# ==============================================================================
# File: 01_data_integration_and_QC.R
# Description: merge data, analyze HLA matrices, perform QC and visualize annotated cells.
# ==============================================================================
source(file.path("D:/项目/2026.Nature Communications/Revison/Code", "Figure4_urine_immune/01_setup_and_helpers.R"))

ensure_dir(FIGURE_DIR)
ensure_dir(RESULT_DIR)
ensure_dir(QC_DIR)
ensure_dir(DEMULTIPLEX_RESULT_DIR)

# Read data and perform initial quality control.
samplenames <- ANALYSIS_BATCHES
scRNA_list <- list()
counter <- 1

# Locate 10x filtered_feature_bc_matrix folders directly under GEM or batch subdirectories.
for (batch in samplenames) {
  tenx_dirs <- find_10x_dirs(batch)
  if (nrow(tenx_dirs) == 0) {
    warning("No 10X matrix found for batch: ", batch)
    next
  }

  for (row_id in seq_len(nrow(tenx_dirs))) {
    data_path <- tenx_dirs$data_dir[row_id]
    run_index <- tenx_dirs$run_index[row_id]
    message("Reading 10X matrix: ", data_path)

    Counts <- Read10X(data.dir = data_path)
    project_name <- paste0(batch, "_", run_index)
    obj <- CreateSeuratObject(counts = Counts, project = project_name, assay = "RNA")
    new_barcodes <- paste0(str_sub(colnames(obj), 1, 16), "-", batch, "-", run_index)
    obj <- RenameCells(obj, new.names = new_barcodes)
    scRNA_list[[counter]] <- obj
    counter <- counter + 1
  }
}
# Merge the Seurat objects.
if (length(scRNA_list) > 1) {
  scRNA_Intergrate <- merge(scRNA_list[[1]], y = scRNA_list[-1], project = "Urine_cell")
  message("Merged ", length(scRNA_list), " samples.")
} else if (length(scRNA_list) == 1) {
  scRNA_Intergrate <- scRNA_list[[1]]
  message("One sample was loaded; merging is not required.")
} else {
  stop("No valid input data found. Check directory structure and paths.")
}

scRNA_Intergrate <- scRNA_QC(scRNA_Intergrate, min_RNA = 200, max_MT = 20, save_path = QC_DIR)

# Identify doublets with scrubletR.
options(future.globals.maxSize = 10000 * 1024^2)

scRNA_Intergrate <- NormalizeData(scRNA_Intergrate) %>% 
  FindVariableFeatures(nfeatures = 2000) %>% 
  ScaleData() %>% RunPCA() %>% 
  RunHarmony(group.by.vars = "orig.ident") %>% 
  JoinLayers()

scRNA_Intergrate <- scrublet(scRNA_Intergrate)
scRNA_Intergrate <- subset(scRNA_Intergrate, subset = doublet_scores < 0.2)

# Export the analysis source data for figure reproduction.
saveRDS(scRNA_Intergrate, object_path("Urine_scRNA_Intergrate.rds"))
