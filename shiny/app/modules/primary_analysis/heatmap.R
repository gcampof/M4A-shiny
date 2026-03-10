suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(ConsensusClusterPlus)
  library(matrixStats)
  library(viridis)
  library(grid)
})

pick_cols_from_vector <- function(levs, base_named_vec,
                                  fallback_seed = 1,
                                  fallback_palette = "Set2") {
  if (!requireNamespace("colorspace", quietly = TRUE)) {
    stop("Package 'colorspace' is required. Install with install.packages('colorspace').")
  }
  
  levs <- as.character(levs)
  out <- base_named_vec[levs]
  missing <- is.na(out)
  
  if (any(missing)) {
    set.seed(fallback_seed)
    out[missing] <- colorspace::qualitative_hcl(sum(missing), palette = fallback_palette)
    names(out)[missing] <- levs[missing]
  }
  
  names(out) <- levs
  out
}

row_cor_dist <- function(mat, method = "pearson") {
  cc <- stats::cor(t(mat), use = "pairwise.complete.obs", method = method)
  as.dist(1 - cc)
}

export_row_names <- function(row_clusters, out_dir, prefix) {
  utils::write.table(
    data.frame(CpG = names(row_clusters),
               row_cluster = as.integer(row_clusters),
               stringsAsFactors = FALSE),
    file = file.path(out_dir, paste0(prefix, "_row_clusters.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}



generate_consensus_clustering <- function ( 
    heatmap_dir,
    beta,
    targets,
    id_col,
    top = 250,
    rowK = 5,
    colK = 7,
    cc_kmax = 9,
    cc_reps = 500,
    cc_pItem = 0.8
){
  
  # Fixed values
  cc_pFeature = 1
  cc_seed = 1
  set.seed(cc_seed)
  prefix_tag = "Panel"
  
  cc_name <- paste0("ConsensusClusterPlus_top", top)
  
  cc <- ConsensusClusterPlus(
    as.matrix(mat),
    maxK = cc_kmax,
    reps = cc_reps,
    pItem = cc_pItem,
    pFeature = cc_pFeature,
    clusterAlg = "hc",
    distance = "pearson",
    innerLinkage = "ward.D2",
    finalLinkage = "ward.D2",
    seed = cc_seed,
    plot = "pdf",
    title = cc_name
  )
  
  col_class <- cc[[colK]]$consensusClass[colnames(mat)]
  col_split <- factor(col_class)
  
  write.table(
    data.frame(sample = names(col_class), CCP_cluster = col_class),
    file = file.path(heatmap_dir, paste0("ConsensusClass_k", colK, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}



prepare_heatmap_data <- function (
    heatmap_dir,
    beta,
    targets,
    id_col,
    sample_group_col,
    annotations = "",
    top_cpgs = 250,
    rowK = 5,
    colK = 7,
    show_row_names = FALSE,
    show_col_names = FALSE){
  
  # Fixed values
  method = "pearson"
  prefix_tag = "Panel"
  pathway_col = "PathwayGene_Affected"
  
  
  # Prepare inputs
  align_res <- align_targets_to_beta_cols(beta, targets, id_col)
  beta2 <- align_res$beta2
  targets2 <- align_res$targets2
  
  # top variable CpGs
  top_cpgs <- min(top_cpgs, nrow(beta2))
  mat <- get_top_mad_probes(beta2, top_cpgs)
  
  top_anno <- HeatmapAnnotation(
    df = data.frame(
      Sample_Group = targets2[[sample_group_col]],
      Pathway = targets2[[pathway_col]]
    ),
    # TODO: maybe remove?
    col = list(
      Sample_Group = sg_cols,
      Pathway = pw_cols
    ),
    annotation_height = unit(c(6,6),"mm")
  )
  
  # Consensus clustering (columns)
  # Maybe here?
  generate_consensus_clustering()

  # =============================================================================
  # Row clustering (1-cor + ward.D2)
  # =============================================================================
  
  row_tree <- hclust(row_cor_dist(mat), method = "ward.D2")
  row_clusters <- cutree(row_tree, k = rowK)
  
  prefix <- paste0(prefix_tag,
                   "_top", top,
                   "_rowK", rowK,
                   "_colK", colK)
  
  # =============================================================================
  # Build Heatmap
  # =============================================================================
  
  ht <- Heatmap(
    mat,
    name = "Methylation",
    col = viridis(100),
    
    cluster_rows = row_tree,
    row_split = rowK,
    
    cluster_columns = TRUE,
    column_split = col_split,
    
    show_row_names = show_row_names,
    show_column_names = show_col_names,
    
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    
    top_annotation = top_anno
  )
  
  pdf(file.path(heat_dir, paste0(prefix, "_heatmap.pdf")), width = 12, height = 8)
  draw(ht, merge_legends = TRUE,
       heatmap_legend_side = "right",
       annotation_legend_side = "bottom")
  dev.off()
  
  invisible(list(
    ht = ht,
    matrix = mat,
    row_clusters = row_clusters,
    col_clusters = col_class,
    heat_dir = heat_dir
  ))
}



