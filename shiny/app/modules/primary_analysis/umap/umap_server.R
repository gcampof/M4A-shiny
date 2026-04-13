source("modules/primary_analysis/umap/umap_utils.R")

umap_server <- function(beta_merged, targets_merged, PALETTES, out_dir, input, output, session) {
  
  umap_data <- reactive({
    req(beta_merged, targets_merged, input$umap_top_cpgs, input$umap_min_dist,
        input$umap_n_neighbors, input$umap_knn, input$umap_consensus_k_max,
        input$umap_metric, input$umap_id_col, input$umap_seed)
    
    validate(need(!is.null(beta_merged), "Beta data missing"))
    validate(need(!is.null(targets_merged), "Targets data missing"))
    
    tryCatch({
      result <- prepare_umap_data(
        beta_merged, 
        targets_merged, 
        input$umap_top_cpgs,
        input$umap_min_dist,
        input$umap_n_neighbors,
        input$umap_metric, 
        input$umap_knn, 
        input$umap_consensus_k_max,
        input$umap_id_col, 
        input$umap_seed
      )

      # Return just the umap_df for plotting
      result$umap_df
      
    }, error = function(e) {
      shiny::validate(
        shiny::need(FALSE, paste0("Error preparing UMAP data: ", e$message))
      )
      NULL
    })
  })
  
  
  output$umap_plot <- renderPlotly({
    # req(umap_data(), input$umap_color_by, input$umap_legend_position,
    #     input$umap_color_palette)
    # 
    validate(need(
      input$umap_color_palette %in% names(PALETTES$all_palettes),
      "Invalid palette"
    ))
    
    tryCatch({
      p <- plot_umap(
        umap_df = umap_data(),
        color_by = input$umap_color_by,
        legend_position = input$umap_legend_position,
        color_palette = PALETTES$all_palettes[[input$umap_color_palette]],
        show_labels = input$umap_show_labels,
        show_summary = input$umap_show_summary,
        top_cpgs = input$umap_top_cpgs,
        min_dist = input$umap_min_dist,
        n_neighbors = input$umap_n_neighbors,
        knn = input$umap_knn,
        consensus_k_max = input$umap_consensus_k_max,
        metric = input$umap_metric,
        out_dir = out_dir
      )
      
      plotly::ggplotly(p, tooltip = c("colour", "label"))
    }, error = function(e) {
      shiny::validate(
        shiny::need(FALSE, paste0("Error rendering UMAP plot: ", e$message))
      )
    })
  })

}