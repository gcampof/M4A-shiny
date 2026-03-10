suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

col_vector <- c(
  "Normal" = "#0d570b",   "JGCT" = "#d97009",  "MST" = "#54c40a",
  "AGCT" = "#52367d",  "Wd-SLCT" = "#bfa21d",  "SLCT" = "#6998b5", "MIC_SRST" ="#E41A1C",
  "HGSC"="#B3E2CD",  "ATRT" = "#e080c3",
  "SCT" = "#874c23",   "SCST" = "#8c8b8b", "SCCOHT"= "#E5C494",
  "Mixed" = "#14b89c",  "LCT" = "#03465e",  "SRST"="#E31A1C", "PDX Granulosa"= "#8c8b8b"
)


align_targets_to_beta_cols <- function(beta, targets, id_col) {
  common_ids <- intersect(colnames(beta), targets[[id_col]])
  if (length(common_ids) < 2) stop("Not enough overlapping sample IDs between beta colnames and targets[[id_col]].")
  
  beta2 <- beta[, common_ids, drop = FALSE]
  targets2 <- targets[match(common_ids, targets[[id_col]]), , drop = FALSE]
  rownames(targets2) <- targets2[[id_col]]
  
  list(beta2 = beta2, targets2 = targets2)
}

prepare_umap_data <- function(
    beta,
    targets,
    top_cpgs,
    min_dist,
    n_neighbors,
    metric,
    id_col
    ) {
  
  # Prepare inputs
  align_res <- align_targets_to_beta_cols(beta, targets, id_col)
  beta2 <- align_res$beta2
  targets2 <- align_res$targets2

  # top variable CpGs
  top_cpgs <- min(top_cpgs, nrow(beta2))
  mat <- get_top_mad_probes(beta2, top_cpgs)
  mat_t <- t(mat)

  # UMAP config
  cfg <- umap::umap.defaults
  cfg$min_dist <- min_dist
  cfg$n_neighbors <- n_neighbors
  cfg$metric <- metric
  
  um <- umap::umap(mat_t, config = cfg)
  
  # output df one row one sample
  df <- data.frame(
    Sample = rownames(um$layout),
    UMAP1 = um$layout[, 1],
    UMAP2 = um$layout[, 2],
    stringsAsFactors = FALSE
  )
  
  # attach metadata
  df <- cbind(df, targets2[df$Sample, , drop = FALSE])

  df
}

plot_umap <- function(
    umap_df,
    color_by,
    show_labels = FALSE,
    show_summary = TRUE,
    top_cpgs,
    min_dist,
    n_neighbors,
    metric
) {
  df <- data.frame(
    X = umap_df$UMAP1,
    Y = umap_df$UMAP2,
    Color = umap_df[[color_by]],
    Label = umap_df$Sample
  )
  
  color_vals <- unique(df$Color)
  matched <- col_vector[names(col_vector) %in% color_vals]
  
  # fallback: generate colors for unmatched values
  unmatched <- setdiff(color_vals, names(col_vector))
  if (length(unmatched) > 0) {
    fallback <- setNames(scales::hue_pal()(length(unmatched)), unmatched)
    matched <- c(matched, fallback)
  }
  
  mode <- if (show_labels) "markers+text" else "markers"
  
  p <- plotly::plot_ly(
    df,
    x = ~X,
    y = ~Y,
    type = "scatter",
    mode = mode,
    color = ~Color,
    colors = matched,
    text = ~Label,
    textposition = "top center",
    marker = list(size = 7, opacity = 0.9)
  )
  
  # ---- summary annotation ----
  if (show_summary) {
    summary_txt <- paste0(
      "CpGs=", top_cpgs,
      " | min_dist=", min_dist,
      " | n_neighbors=", n_neighbors,
      " | metric=", metric,
      " | color=", color_by,
      if (show_labels) " (labeled)" else ""
    )
    
    p <- p %>% plotly::layout(
      annotations = list(
        list(
          x = 0,
          y = 2,
          xref = "paper",
          yref = "paper",
          xanchor = "left",
          yanchor = "top",
          showarrow = FALSE,
          text = summary_txt,
          font = list(size = 12)
        )
      )
    )
  }
  
  p %>% plotly::layout(
    xaxis = list(title = "UMAP1"),
    yaxis = list(title = "UMAP2")
  )
}
