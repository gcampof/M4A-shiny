# Get top cpg
get_top_mad_probes <- function(beta, n) {
  beta <- as.matrix(beta)
  storage.mode(beta) <- "numeric"
  
  idx <- order(matrixStats::rowMads(beta), decreasing = TRUE)[seq_len(n)]
  beta[idx, , drop = FALSE]
}


# Prepare MDS data from beta matrix
prepare_mds_data <- function(beta_merged, top_cpgs) {
  # beta_merged: first column = CpG ID, rest = samples
  message("calculating mds")
  beta <- as.matrix(beta_merged[, -1])
  rownames(beta) <- beta_merged[[1]]
  storage.mode(beta) <- "numeric"
  
  # 1. select top variable CpGs
  beta_top <- get_top_mad_probes(beta, top_cpgs)

  # 2. transpose: samples x CpGs
  beta_t <- t(beta_top)
  
  # 3. remove zero-variance samples (important!)
  sds <- matrixStats::rowSds(beta_t)
  beta_t <- beta_t[sds > 0, , drop = FALSE]
  
  req(nrow(beta_t) >= 3)
  
  # 4. distance + MDS
  dist_mat <- dist(beta_t, method = "euclidean")
  coords <- cmdscale(dist_mat, k = 2)

  # 5. return clean data.frame
  data.frame(
    Sample = rownames(coords),
    Dim1 = coords[, 1],
    Dim2 = coords[, 2],
    stringsAsFactors = FALSE
  )
}


# Create MDS plotly plot
plot_mds <- function(mds_df, meta, color_by) {

  plotly::plot_ly(
    data = mds_df,
    x = ~Dim1,
    y = ~Dim2,
    type = "scatter",
    mode = "markers",
    color = meta[[color_by]],
    text = ~Sample
  ) %>%
    plotly::layout(
      xaxis = list(title = "MDS 1"),
      yaxis = list(title = "MDS 2")
    )
}



# Prepare PCA data from beta matrix
prepare_pca_data <- function(beta_merged, top_cpgs) {
  message("calculating pca")
  
  beta <- as.matrix(beta_merged[, -1])
  rownames(beta) <- beta_merged[[1]]
  storage.mode(beta) <- "numeric"
  
  beta_top <- get_top_mad_probes(beta, top_cpgs)
  
  # remove zero-variance CpGs
  sds_cpg <- matrixStats::rowSds(beta_top)
  beta_top <- beta_top[sds_cpg > 0, , drop = FALSE]
  
  # remove zero-variance samples
  sds_samples <- matrixStats::colSds(beta_top)
  beta_top <- beta_top[, sds_samples > 0, drop = FALSE]
  
  req(ncol(beta_top) >= 3)
  
  pca <- prcomp(t(beta_top), scale. = TRUE)
 
  pca
}


# Create PCA plotly plot
plot_pca <- function(pca, meta, color_by, dims = "PC1 vs PC2") {

  dim_indices <- switch(
    dims,
    "PC1 vs PC2" = c(1, 2),
    "PC1 vs PC3" = c(1, 3),
    "PC2 vs PC3" = c(2, 3)
  )
  
  pcs <- pca$x
  
  df <- data.frame(
    Sample = rownames(pcs),
    X = pcs[, dim_indices[1]],
    Y = pcs[, dim_indices[2]],
    Color = meta[[color_by]]
  )
  
  plotly::plot_ly(
    df,
    x = ~X,
    y = ~Y,
    type = "scatter",
    mode = "markers",
    color = ~Color,
    text = ~Sample
  ) %>%
    plotly::layout(
      xaxis = list(title = paste0("PC", dim_indices[1])),
      yaxis = list(title = paste0("PC", dim_indices[2]))
    )
}