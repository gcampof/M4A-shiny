library(plotly)
library(RColorBrewer)

generate_beta_boxplots <- function(array, beta, results_path) {
  message("Calculating boxplot statistics...")
  stats_list <- lapply(colnames(beta), function(sample) {
    values <- beta[, sample]
    q1 <- quantile(values, 0.25, na.rm = TRUE)
    q2 <- quantile(values, 0.5, na.rm = TRUE)
    q3 <- quantile(values, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    
    lower_limit <- q1 - 1.5 * iqr
    upper_limit <- q3 + 1.5 * iqr
    
    whisker_low <- min(values[values >= lower_limit], na.rm = TRUE)
    whisker_high <- max(values[values <= upper_limit], na.rm = TRUE)
    
    data.frame(
      sample = sample,
      q1 = as.numeric(q1),
      median = as.numeric(q2),
      q3 = as.numeric(q3),
      whisker_low = as.numeric(whisker_low),
      whisker_high = as.numeric(whisker_high),
      stringsAsFactors = FALSE
    )
  })
  
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL
  
  # Generar paleta de colores
  n_samples <- nrow(stats_df)
  pal <- RColorBrewer::brewer.pal(min(n_samples, 8), "Set2")
  colors <- rep(pal, length.out = n_samples)
  
  p <- plot_ly(stats_df, x = ~factor(sample, levels = sample)) %>%
    add_boxplot(
      q1 = ~q1,
      median = ~median,
      q3 = ~q3,
      lowerfence = ~whisker_low,
      upperfence = ~whisker_high,
      fillcolor = ~factor(sample, levels = sample),
      line = list(color = "gray50"),
      marker = list(opacity = 0),
      boxmean = FALSE,
      hovertemplate = "<b>%{x}</b><br>Q1: %{q1:.4f}<br>Median: %{median:.4f}<br>Q3: %{q3:.4f}<extra></extra>",
      showlegend = FALSE
    ) %>%
    layout(
      title = list(
        text = paste("Beta value distribution –", array),
        font = list(size = 16)
      ),
      xaxis = list(
        title = "Sample",
        tickangle = -45,
        automargin = TRUE
      ),
      yaxis = list(
        title = "Beta value",
        gridcolor = "rgb(240, 240, 240)"
      ),
      template = "plotly_white",
      hovermode = "closest",
      margin = list(b = 100, l = 60, r = 40, t = 80),
      height = 600,
      width = 1000,
      font = list(size = 12, family = "sans-serif"),
      plot_bgcolor = "rgba(255, 255, 255, 1)",
      paper_bgcolor = "rgba(255, 255, 255, 1)"
    )
  
  # Guardar HTML
  results_dir <- file.path(results_path, array)
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  html_path <- file.path(results_dir, paste0("beta_boxplot_", array, ".html"))
  htmlwidgets::saveWidget(p, html_path, selfcontained = TRUE)
  
  return(p)
}