# Purpose: Upstream BALF processing, repertoire and disease-signature computation.
# Original source: run_sc_pseudobulk_functional_analysis.R
# Data cohort: balf. Inputs use the compact English data layout.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(DESeq2)
  library(clusterProfiler)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
})

# ==============================================================================
# Module 1: infrastructure and preprocessing.
# ==============================================================================

# Load the required species annotation database.
.get_orgdb <- function(organism = c("human", "mouse")) {
  organism <- match.arg(organism)
  
  if (organism == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Install org.Hs.eg.db first: BiocManager::install('org.Hs.eg.db')")
    }
    OrgDb <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    kegg_organism <- "hsa"
  } else if (organism == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Install org.Mm.eg.db first: BiocManager::install('org.Mm.eg.db')")
    }
    OrgDb <- getExportedValue("org.Mm.eg.db", "org.Mm.eg.db")
    kegg_organism <- "mmu"
  }
  
  list(OrgDb = OrgDb, kegg_organism = kegg_organism)
}

# Retrieve counts compatibly with Seurat v4/v5.
.get_counts <- function(seu, assay = "RNA", layer = "counts") {
  tryCatch(
    { SeuratObject::GetAssayData(seu, assay = assay, layer = layer) },
    error = function(e) { SeuratObject::GetAssayData(seu, assay = assay, slot = layer) }
  )
}

# Construct pseudobulk count matrices and metadata.
make_pseudobulk_counts <- function(
    seu,
    group_col = "Group",
    sample_col = "HLA_Splitter",
    celltype_col = NULL,
    default_celltype = NULL,
    assay = "RNA",
    layer = "counts"
) {
  counts <- .get_counts(seu, assay = assay, layer = layer)
  meta <- seu@meta.data[colnames(counts), , drop = FALSE]
  
  # Check required metadata columns.
  required_cols <- c(group_col, sample_col)
  if (!is.null(celltype_col)) required_cols <- c(required_cols, celltype_col)
  
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Standardize metadata column names and fill configured defaults.
  meta$.group  <- as.character(meta[[group_col]])
  meta$.sample <- as.character(meta[[sample_col]])
  meta$.celltype <- if (is.null(celltype_col)) default_celltype else as.character(meta[[celltype_col]])
  
  # Exclude cells with invalid metadata.
  keep_cells <- complete.cases(meta[, c(".group", ".sample", ".celltype")])
  meta <- meta[keep_cells, , drop = FALSE]
  counts <- counts[, rownames(meta), drop = FALSE]
  
  # Check that each sample maps to one clinical group.
  sample_group_check <- meta %>%
    dplyr::distinct(.sample, .group) %>%
    dplyr::count(.sample) %>%
    dplyr::filter(n > 1)
  
  if (nrow(sample_group_check) > 0) {
    stop("Samples map to multiple Group values; check metadata. Affected samples: ", 
         paste(sample_group_check$.sample, collapse = ", "))
  }
  
  # Aggregate counts by sparse-matrix multiplication.
  meta$.pb_id <- paste(meta$.celltype, meta$.sample, sep = "__")
  pb_factor <- factor(meta$.pb_id)
  
  design_mat <- Matrix::sparse.model.matrix(~ 0 + pb_factor)
  colnames(design_mat) <- levels(pb_factor)
  
  pb_counts <- as(counts %*% design_mat, "dgCMatrix")
  
  pb_meta <- meta %>%
    dplyr::group_by(.pb_id) %>%
    dplyr::summarise(
      celltype = dplyr::first(.celltype),
      sample   = dplyr::first(.sample),
      group    = dplyr::first(.group),
      n_cells  = dplyr::n(),
      .groups  = "drop"
    ) %>%
    as.data.frame() %>%
    tibble::column_to_rownames(".pb_id")
  
  # Align pseudobulk matrix columns with metadata rows.
  pb_counts <- pb_counts[, rownames(pb_meta), drop = FALSE]
  
  list(counts = pb_counts, metadata = pb_meta)
}


# ==============================================================================
# Module 2: core statistics and functional annotation.
# ==============================================================================

# Define the gene-identifier conversion helper.
map_genes_to_entrez <- function(genes, OrgDb, from_type = "SYMBOL") {
  genes <- unique(genes[complete.cases(genes) & genes != ""])
  if (length(genes) == 0) return(data.frame())
  
  suppressMessages({
    mp <- clusterProfiler::bitr(
      genes, fromType = from_type, toType = "ENTREZID", OrgDb = OrgDb
    )
  })
  
  if (is.null(mp) || nrow(mp) == 0) return(data.frame())
  unique(mp[, c(from_type, "ENTREZID")])
}

# Run per-cell-type DESeq2 analysis with global LRT and pairwise Wald tests.
run_deseq_for_one_celltype <- function(
    pb_counts, pb_meta, celltype, comparisons,
    covariates = NULL, min_cells_per_sample = 20,
    min_samples_per_group = 2, min_count = 10, min_samples_expr = 2
) {
  # Filter samples.
  meta_ct <- pb_meta %>% filter(celltype == !!celltype, n_cells >= min_cells_per_sample)
  if (nrow(meta_ct) == 0) {
    message("Skipping ", celltype, ": insufficient samples")
    return(NULL)
  }
  
  group_tab <- table(meta_ct$group)
  if (sum(group_tab >= min_samples_per_group) < 2) {
    message("Skipping ", celltype, ": fewer than two eligible groups")
    return(NULL)
  }
  
  # Prepare data and filter genes.
  cts <- as.matrix(pb_counts[, rownames(meta_ct), drop = FALSE])
  storage.mode(cts) <- "integer"
  meta_ct$group <- factor(meta_ct$group)
  
  keep_gene <- rowSums(cts >= min_count) >= min_samples_expr
  cts <- cts[keep_gene, , drop = FALSE]
  if (nrow(cts) < 100) {
    message("Skipping ", celltype, ": too few genes remain after filtering")
    return(NULL)
  }
  
  # Construct the model formula.
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, colnames(meta_ct))
    if (length(missing_cov) > 0) stop("Missing covariates: ", paste(missing_cov, collapse = ", "))
  }
  
  full_formula <- as.formula(paste("~", paste(c(covariates, "group"), collapse = " + ")))
  reduced_formula <- as.formula(if(length(covariates) > 0) paste("~", paste(covariates, collapse = " + ")) else "~ 1")
  
  # Run the likelihood-ratio test.
  overall_res <- tryCatch({
    dds_lrt <- DESeq2::DESeqDataSetFromMatrix(cts, colData = meta_ct, design = full_formula)
    dds_lrt <- DESeq2::DESeq(dds_lrt, test = "LRT", reduced = reduced_formula, quiet = TRUE)
    as.data.frame(DESeq2::results(dds_lrt)) %>%
      rownames_to_column("gene") %>%
      as_tibble() %>%
      mutate(celltype = celltype, test = "LRT_group_overall") %>%
      relocate(celltype, test, gene)
  }, error = function(e) { message("LRT failed: ", e$message); tibble() })
  
  # Run Wald tests.
  dds_wald <- tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(cts, colData = meta_ct, design = full_formula)
    DESeq2::DESeq(dds, quiet = TRUE)
  }, error = function(e) { message("DESeq2 Wald test failed: ", e$message); NULL })
  
  pairwise_list <- list()
  if (!is.null(dds_wald)) {
    for (comp in comparisons) {
      num <- comp[1]; den <- comp[2]
      contrast_name <- paste0(num, "_vs_", den)
      
      sub_group <- as.character(meta_ct$group[meta_ct$group %in% c(num, den)])
      
      if (!all(c(num, den) %in% sub_group) || any(table(sub_group) < min_samples_per_group)) {
        message("Skipping ", contrast_name, ": insufficient eligible samples for this comparison.")
        next
      }
      
      res <- tryCatch({
        as.data.frame(DESeq2::results(dds_wald, contrast = c("group", num, den))) %>%
          rownames_to_column("gene") %>%
          as_tibble() %>%
          mutate(celltype = celltype, contrast = contrast_name, numerator = num, denominator = den) %>%
          relocate(celltype, contrast, numerator, denominator, gene) %>%
          arrange(padj)
      }, error = function(e) { tibble() })
      
      pairwise_list[[contrast_name]] <- res
    }
  }
  
  list(overall = overall_res, pairwise = pairwise_list, dds = dds_wald, metadata = meta_ct)
}

# Run GO/KEGG over-representation and GSEA for one contrast.
run_functional_for_contrast <- function(
    deg_table, OrgDb, kegg_organism = "hsa", gene_id_type = "SYMBOL",
    padj_cutoff = 0.05, lfc_cutoff = 0.25, go_ont = "BP",
    min_gs_size = 10, max_gs_size = 500
) {
  deg_table <- deg_table %>% filter(!is.na(gene))
  
  # Extract background genes and differential genes.
  universe_entrez <- unique(map_genes_to_entrez(deg_table$gene[!is.na(deg_table$pvalue)], OrgDb, gene_id_type)$ENTREZID)
  
  get_entrez <- function(direction_condition) {
    genes <- deg_table %>% filter(!!rlang::parse_expr(direction_condition)) %>% pull(gene) %>% unique()
    unique(map_genes_to_entrez(genes, OrgDb, gene_id_type)$ENTREZID)
  }
  
  up_entrez   <- get_entrez(paste0("!is.na(padj) & padj < ", padj_cutoff, " & log2FoldChange > ", lfc_cutoff))
  down_entrez <- get_entrez(paste0("!is.na(padj) & padj < ", padj_cutoff, " & log2FoldChange < -", lfc_cutoff))
  
  safe_call <- function(expr) tryCatch(expr, error = function(e) NULL)
  
  # Run over-representation analysis.
  run_ora <- function(genes, func, ...) {
    if (length(genes) >= 5) safe_call(func(gene = genes, universe = universe_entrez, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, ...)) else NULL
  }
  
  ego_up   <- run_ora(up_entrez, clusterProfiler::enrichGO, OrgDb = OrgDb, keyType = "ENTREZID", ont = go_ont, readable = TRUE)
  ego_down <- run_ora(down_entrez, clusterProfiler::enrichGO, OrgDb = OrgDb, keyType = "ENTREZID", ont = go_ont, readable = TRUE)
  ekegg_up   <- run_ora(up_entrez, clusterProfiler::enrichKEGG, organism = kegg_organism)
  ekegg_down <- run_ora(down_entrez, clusterProfiler::enrichKEGG, organism = kegg_organism)
  
  # Run gene-set enrichment analysis.
  rank_df <- deg_table %>%
    mutate(score = case_when(
      !is.na(stat) ~ stat,
      is.na(stat) & !is.na(pvalue) & !is.na(log2FoldChange) ~ sign(log2FoldChange) * (-log10(pmax(pvalue, .Machine$double.xmin))),
      TRUE ~ NA_real_
    )) %>%
    filter(!is.na(score)) %>% select(gene, score)
  
  rank_map <- map_genes_to_entrez(rank_df$gene, OrgDb, gene_id_type)
  
  gene_list <- rank_df %>%
    inner_join(rank_map, by = c("gene" = gene_id_type)) %>%
    arrange(desc(abs(score))) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    { setNames(.$score, .$ENTREZID) } %>%
    sort(decreasing = TRUE)
  
  gsea_go <- gsea_kegg <- NULL
  if (length(gene_list) >= 100) {
    gsea_go   <- safe_call(clusterProfiler::gseGO(geneList = gene_list, OrgDb = OrgDb, keyType = "ENTREZID", ont = go_ont, minGSSize = min_gs_size, maxGSSize = max_gs_size, pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE))
    gsea_kegg <- safe_call(clusterProfiler::gseKEGG(geneList = gene_list, organism = kegg_organism, minGSSize = min_gs_size, maxGSSize = max_gs_size, pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE))
  }
  
  list(
    GO_UP = ego_up, GO_DOWN = ego_down, KEGG_UP = ekegg_up, KEGG_DOWN = ekegg_down,
    GSEA_GO = gsea_go, GSEA_KEGG = gsea_kegg
  )
}

# Convert enrichment results to a long data frame.
.enrich_to_table <- function(enrich_obj, method, direction, celltype, contrast) {
  if (is.null(enrich_obj)) return(tibble())
  df <- tryCatch(as_tibble(as.data.frame(enrich_obj)), error = function(e) tibble())
  if (nrow(df) == 0) return(tibble())
  
  df %>%
    mutate(
      direction = if ("NES" %in% colnames(.)) ifelse(NES > 0, "UP_IN_NUMERATOR", "DOWN_IN_NUMERATOR") else direction,
      celltype = celltype, contrast = contrast, method = method
    ) %>%
    relocate(celltype, contrast, method, direction)
}


# ==============================================================================
# Module 3: downstream pattern analysis.
# ==============================================================================

# Extract the pathway heatmap matrix.
make_pathway_direction_matrix <- function(pathway_table, method_keep = c("GSEA_GO", "GSEA_KEGG"), padj_cutoff = 0.05) {
  pathway_table %>%
    filter(method %in% method_keep, !is.na(p.adjust), p.adjust < padj_cutoff, !is.na(NES)) %>%
    mutate(
      pathway = paste(method, ID, Description, sep = "__"),
      celltype_contrast = paste(celltype, contrast, sep = "__")
    ) %>%
    select(pathway, celltype_contrast, NES) %>%
    distinct() %>%
    pivot_wider(names_from = celltype_contrast, values_from = NES, values_fill = 0)
}

# Identify pathways with reversed direction across comparisons.
find_reversed_pathways <- function(pathway_table, contrast_a, contrast_b, method_keep = c("GSEA_GO", "GSEA_KEGG"), padj_cutoff = 0.05) {
  extract_contrast <- function(ctr) {
    pathway_table %>%
      filter(contrast == ctr, method %in% method_keep, !is.na(p.adjust), p.adjust < padj_cutoff, !is.na(NES)) %>%
      select(celltype, method, ID, Description, NES, padj = p.adjust)
  }
  
  inner_join(
    extract_contrast(contrast_a) %>% rename(NES_a = NES, padj_a = padj),
    extract_contrast(contrast_b) %>% rename(NES_b = NES, padj_b = padj),
    by = c("celltype", "method", "ID", "Description")
  ) %>%
    filter(sign(NES_a) != sign(NES_b)) %>%
    mutate(reversal_score = abs(NES_a) + abs(NES_b)) %>%
    arrange(desc(reversal_score))
}

# Identify cell-type-specific pathways.
find_celltype_specific_pathways <- function(pathway_table, contrast_name, method_keep = c("GSEA_GO", "GSEA_KEGG", "GO_ORA", "KEGG_ORA"), padj_cutoff = 0.05) {
  pathway_table %>%
    filter(contrast == contrast_name, method %in% method_keep, !is.na(p.adjust), p.adjust < padj_cutoff) %>%
    group_by(method, ID, Description) %>%
    mutate(n_celltypes = n_distinct(celltype)) %>%
    ungroup() %>%
    filter(n_celltypes == 1) %>%
    arrange(celltype, method, p.adjust)
}


# ==============================================================================
# Module 4: main analysis wrapper.
# ==============================================================================

# Define the main analysis function.
run_sc_pseudobulk_functional_analysis <- function(
    seu, group_col = "Group", sample_col = "HLA_Splitter", celltype_col = NULL, default_celltype = NULL,
    comparisons = NULL, organism = c("human", "mouse"),
    gene_id_type = "SYMBOL", assay = "RNA", layer = "counts", covariates = NULL,
    min_cells_per_sample = 20, min_samples_per_group = 2, min_count = 10, min_samples_expr = 2,
    padj_cutoff = 0.05, lfc_cutoff = 0.25, go_ont = "BP", output_dir = NULL
) {
  # Construct the analysis configuration.
  org_info <- .get_orgdb(organism)
  pb <- make_pseudobulk_counts(seu, group_col, sample_col, celltype_col, default_celltype, assay, layer)
  
  celltypes <- unique(pb$metadata$celltype)
  all_overall <- list(); all_pairwise <- list(); functional_objects <- list(); functional_tables <- list()
  
  # Run the analysis for each cell type.
  for (ct in celltypes) {
    message("========== Analyzing cell type: ", ct, " ==========")
    de_ct <- run_deseq_for_one_celltype(pb$counts, pb$metadata, ct, comparisons, covariates, min_cells_per_sample, min_samples_per_group, min_count, min_samples_expr)
    if (is.null(de_ct)) next
    
    all_overall[[ct]] <- de_ct$overall
    functional_objects[[ct]] <- list()
    
    for (contrast_name in names(de_ct$pairwise)) {
      deg_tbl <- de_ct$pairwise[[contrast_name]]
      all_pairwise[[paste(ct, contrast_name, sep = "__")]] <- deg_tbl
      
      if (nrow(deg_tbl) > 0) {
        message("  -> Functional enrichment: ", contrast_name)
        fun_res <- run_functional_for_contrast(deg_tbl, org_info$OrgDb, org_info$kegg_organism, gene_id_type, padj_cutoff, lfc_cutoff, go_ont)
        
        functional_objects[[ct]][[contrast_name]] <- fun_res
        functional_tables[[paste(ct, contrast_name, sep = "__")]] <- bind_rows(
          .enrich_to_table(fun_res$GO_UP, "GO_ORA", "UP", ct, contrast_name),
          .enrich_to_table(fun_res$GO_DOWN, "GO_ORA", "DOWN", ct, contrast_name),
          .enrich_to_table(fun_res$KEGG_UP, "KEGG_ORA", "UP", ct, contrast_name),
          .enrich_to_table(fun_res$KEGG_DOWN, "KEGG_ORA", "DOWN", ct, contrast_name),
          .enrich_to_table(fun_res$GSEA_GO, "GSEA_GO", "GSEA", ct, contrast_name),
          .enrich_to_table(fun_res$GSEA_KEGG, "GSEA_KEGG", "GSEA", ct, contrast_name)
        )
      }
    }
  }
  
  # Summarize cell proportions.
  cell_proportion <- NULL
  if (!is.null(celltype_col)) {
    cell_proportion <- seu@meta.data %>%
      dplyr::filter(!is.na(!!sym(group_col)), !is.na(!!sym(sample_col)), !is.na(!!sym(celltype_col))) %>%
      dplyr::count(sample = !!sym(sample_col), group = !!sym(group_col), celltype = !!sym(celltype_col), name = "n_cells") %>%
      dplyr::group_by(sample, group) %>%
      dplyr::mutate(total_cells = sum(n_cells), proportion = n_cells / total_cells) %>%
      dplyr::ungroup()
  }
  
  # Combine and save result tables.
  overall_table <- bind_rows(all_overall)
  pairwise_table <- bind_rows(all_pairwise)
  pathway_table <- bind_rows(functional_tables)
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(overall_table, file.path(output_dir, "01_overall_LRT_group_test.csv"))
    readr::write_csv(pairwise_table, file.path(output_dir, "02_pairwise_DEG_all_contrasts.csv"))
    readr::write_csv(pathway_table, file.path(output_dir, "03_GO_KEGG_GSEA_pathway_table.csv"))
    if (!is.null(cell_proportion)) readr::write_csv(cell_proportion, file.path(output_dir, "04_celltype_proportion_by_sample.csv"))
  }
  
  list(
    pseudobulk = pb, overall_LRT = overall_table, pairwise_DEG = pairwise_table,
    functional_objects = functional_objects, pathway_table = pathway_table, cell_proportion = cell_proportion
  )
}

# ==============================================================================
# Module 5: cluster dynamic-gene patterns and draw LRT heatmaps/trends.
# ==============================================================================

# Internal helper: retrieve a cell-type-specific z-score matrix.
.get_scaled_lrt_matrix <- function(
    lrt_table, pb_counts, pb_meta, target_celltype, 
    filter_by = c("padj", "pvalue"), p_cutoff = 0.01, top_n = 1000
) {
  # Validate and match argument values.
  filter_by <- match.arg(filter_by)
  
  # Select columns used for filtering and ranking.
  sig_genes <- lrt_table %>%
    dplyr::filter(
      celltype == target_celltype, 
      !is.na(.data[[filter_by]]), 
      .data[[filter_by]] < p_cutoff
    ) %>%
    dplyr::arrange(.data[[filter_by]]) %>%
    head(top_n) %>%
    dplyr::pull(gene)
  
  if (length(sig_genes) < 10) {
    stop(
      target_celltype, " has fewer than 10 genes satisfying ", filter_by, " < ", p_cutoff, 
      ". Plot skipped.\nReview p_cutoff or filter_by before changing the analysis criteria."
    )
  }
  
  # Select samples belonging to the requested cell type.
  meta_ct <- pb_meta %>%
    dplyr::filter(celltype == target_celltype) %>%
    dplyr::arrange(group) 
  
  # Extract counts and transform to log2(CPM + 1).
  cts <- pb_counts[sig_genes, rownames(meta_ct), drop = FALSE]
  cpm_mat <- sweep(cts, 2, colSums(cts), FUN = "/") * 1e6
  logcpm_mat <- log2(cpm_mat + 1)
  
  # Scale rows to z scores to compare relative expression patterns.
  z_mat <- t(scale(t(logcpm_mat)))
  z_mat <- z_mat[complete.cases(z_mat), , drop = FALSE]
  
  list(
    z_matrix = z_mat,
    metadata = meta_ct,
    filter_metric = filter_by
  )
}


# Draw the global LRT heatmap.
plot_lrt_global_heatmap <- function(
    lrt_table, pb_counts, pb_meta, target_celltype, 
    filter_by = "padj", p_cutoff = 0.01, top_n = 1000, 
    k_clusters = 4, show_rownames = FALSE
) {
  # Retrieve the z-score matrix.
  mat_info <- .get_scaled_lrt_matrix(lrt_table, pb_counts, pb_meta, target_celltype, filter_by, p_cutoff, top_n)
  z_mat <- mat_info$z_matrix
  meta_ct <- mat_info$metadata
  metric <- mat_info$filter_metric
  
  # Prepare annotation tracks.
  annotation_col <- data.frame(
    Group = meta_ct$group,
    row.names = rownames(meta_ct)
  )
  
  # Draw the heatmap with k-means row splits.
  p <- pheatmap::pheatmap(
    mat = z_mat,
    cluster_cols = FALSE,     
    cluster_rows = TRUE,      
    kmeans_k = NA,            
    cutree_rows = k_clusters, 
    annotation_col = annotation_col,
    show_colnames = FALSE,
    show_rownames = show_rownames,
    scale = "none",           
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    main = paste0(target_celltype, " - Dynamic Genes (", metric, " < ", p_cutoff, ")\nTop ", nrow(z_mat), " genes")
  )
  
  return(p)
}


# Extract gene modules and plot their expression trends.
plot_lrt_cluster_trends <- function(
    lrt_table, pb_counts, pb_meta, target_celltype, 
    filter_by = "padj", p_cutoff = 0.01, top_n = 1000, 
    k_clusters = 4
) {
  mat_info <- .get_scaled_lrt_matrix(lrt_table, pb_counts, pb_meta, target_celltype, filter_by, p_cutoff, top_n)
  z_mat <- mat_info$z_matrix
  meta_ct <- mat_info$metadata
  metric <- mat_info$filter_metric
  
  # Cluster gene expression patterns using k-means.
  set.seed(42) 
  km <- kmeans(z_mat, centers = k_clusters, iter.max = 50)
  
  # Reshape z scores to long format and join clinical groups.
  long_df <- as.data.frame(z_mat) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(cluster = paste0("Cluster ", km$cluster[gene])) %>%
    tidyr::pivot_longer(
      cols = -c(gene, cluster), 
      names_to = ".pb_id",    # Name the reshaped sample identifier column .pb_id.
      values_to = "z_score"
    ) %>%
    dplyr::left_join(
      meta_ct %>% tibble::rownames_to_column(".pb_id") %>% dplyr::select(.pb_id, group), 
      by = ".pb_id"
    )
  
  # Calculate group means and confidence intervals for each cluster.
  trend_df <- long_df %>%
    dplyr::group_by(cluster, group) %>%
    dplyr::summarise(
      mean_z = mean(z_score),
      sd_z = sd(z_score),
      n = dplyr::n(),
      se_z = sd_z / sqrt(n),
      .groups = "drop"
    ) %>%
    dplyr::mutate(group = factor(group)) 
  
  # Plot expression trends using ggplot2.
  p <- ggplot(trend_df, aes(x = group, y = mean_z, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(color = "#0073C2FF", linewidth = 1.2) +
    geom_point(color = "#0073C2FF", size = 3) +
    geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z), width = 0.2, color = "#0073C2FF") +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(target_celltype, " - Gene Expression Patterns"),
      subtitle = paste0("K-means clustering (k=", k_clusters, ") | Filter: ", metric, " < ", p_cutoff),
      x = "Disease Group",
      y = "Mean Z-score of Expression"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(face = "bold")
    )
  
  gene_clusters <- data.frame(
    gene = names(km$cluster),
    cluster = paste0("Cluster ", km$cluster),
    row.names = NULL
  )
  
  list(
    plot = p,
    gene_clusters = gene_clusters
  )
}
