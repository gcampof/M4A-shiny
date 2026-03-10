plot_global_methylation <- function(
    beta,
    targets,
    id_col,
    comparison_col,
    comparison_type = "between", 
    group1          = NULL,      
    group2          = NULL,   
    annot,
    color_palette
) {
  # Prepare inputs
  align_res <- align_targets_to_beta_cols(beta, targets, id_col)
  beta2    <- align_res$beta2
  targets2 <- align_res$targets2
  
  # Extract and clean groups
  groups <- trimws(as.character(targets2[[comparison_col]]))
  groups[groups == ""] <- NA
  
  #NA filtering
  keep     <- !is.na(groups)
  targets2 <- targets2[keep, , drop = FALSE]
  beta2    <- beta2[, keep, drop = FALSE]
  groups   <- groups[keep]
  
  # Handle comparison type
  if (comparison_type == "custom") {
    
    shiny::validate(
      shiny::need(length(group1) > 0, "Please assign at least one level to Group 1."),
      shiny::need(length(group2) > 0, "Please assign at least one level to Group 2."),
      shiny::need(
        length(intersect(group1, group2)) == 0,
        paste0("Levels cannot be in both groups: ",
               paste(intersect(group1, group2), collapse = ", "))
      )
    )
    
    # Recode groups: levels in group1 -> "Group 1", group2 -> "Group 2", rest dropped
    keep2  <- groups %in% c(group1, group2)
    beta2  <- beta2[, keep2, drop = FALSE]
    groups <- ifelse(groups[keep2] %in% group1, "Group 1", "Group 2")
    
    comparison_label <- paste0(
      "Group 1: [", paste(group1, collapse = ", "), "]  vs  ",
      "Group 2: [", paste(group2, collapse = ", "), "]"
    )
    comparisons_list <- list(c("Group 1", "Group 2"))
    
  } else {
    
    shiny::validate(
      shiny::need(length(unique(groups)) >= 2,
                  paste0("Column '", comparison_col, "' needs at least 2 groups."))
    )
    
    comparison_label <- comparison_col
    present_groups   <- sort(unique(groups))
    comparisons_list <- combn(present_groups, 2, simplify = FALSE)
  }
  
  # Annotation subsets
  beta_all    <- beta2
  beta_cgi    <- beta2[rownames(beta2) %in% rownames(annot[annot$Relation_to_Island == "Island", ]),]
  beta_shore  <- beta2[rownames(beta2) %in% rownames(annot[annot$Relation_to_Island %in% c("N_Shore", "S_Shore"), ]),]
  beta_shelves <- beta2[rownames(beta2) %in% rownames(annot[annot$Relation_to_Island %in% c("S_Shelf", "N_Shelf", "OpenSea"), ]),]
  
  get_means <- function(mat, label) {
    data.frame(
      SampleMean = colMeans(mat, na.rm = TRUE),
      Group      = groups, 
      CpGType    = label,
      stringsAsFactors = FALSE
    )
  }
  
  plot_data <- rbind(
    get_means(beta_all, "All CpGs"),
    get_means(beta_cgi, "CGI"),
    get_means(beta_shore, "Shores"),
    get_means(beta_shelves, "Shelves & OpenSea")
  )
  
  plot_data$CpGType <- factor(plot_data$CpGType,
                              levels = c("All CpGs", "CGI", "Shores", "Shelves & OpenSea"))
  plot_data$Group   <- factor(plot_data$Group)
  
  # Color mapping
  group_levels <- levels(plot_data$Group)
  n_needed     <- length(group_levels)
  base_colors  <- color_palette(n_needed)
  if (length(base_colors) < n_needed) base_colors <- rep_len(base_colors, n_needed)
  fill_colors  <- setNames(base_colors, group_levels)
  
  # Plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(Group, SampleMean, fill = Group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, linewidth = 0.5) +
    ggplot2::geom_jitter(width = 0.2, size = 1.2, alpha = 0.7) +
    ggpubr::stat_compare_means(
      comparisons = comparisons_list,
      method      = "wilcox.test",
      label       = "p.signif",
      size        = 2.5
    ) +
    ggplot2::facet_wrap(~CpGType, nrow = 1) +
    ggplot2::scale_fill_manual(values = fill_colors) +
    ggplot2::coord_cartesian(ylim = c(0.18, 0.80)) +
    ggplot2::labs(
      title = paste("Global Mean Methylation —", comparison_label),
      x     = NULL,
      y     = "Mean beta per sample"
    ) +
    ggplot2::theme_gray(base_size = 12) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 13, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 11),
      axis.text.y  = ggplot2::element_text(size = 10),
      legend.text      = ggplot2::element_text(size = 12),
      legend.title     = ggplot2::element_text(size = 13, face = "bold"),
      legend.key.size  = ggplot2::unit(1.2, "cm"),    
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid   = ggplot2::element_blank(),
      strip.text   = ggplot2::element_text(size = 11, face = "bold")  # facet labels
    )
  
  p
}