source("modules/primary_analysis/pca/pca_utils.R")
pca_server <- function(beta_merged, targets_merged, PALETTES, out_dir, input, output, session) {
  
  pca_data <- reactive({
    req(beta_merged, targets_merged, input$pca_id_col, input$pca_top_cpgs)
    
    validate(need(!is.null(beta_merged), "Beta data missing"))
    validate(need(!is.null(targets_merged), "Targets data missing"))
    
    tryCatch({
      prepare_pca_data(beta_merged, targets_merged, input$pca_id_col, input$pca_top_cpgs)
    }, error = function(e) {
      shiny::validate(
        shiny::need(FALSE, paste0("Error preparing PCA data: ", e$message))
      )
      NULL
    })
  })
  
  output$pca_plot <- renderPlotly({
    req(pca_data(), input$pca_color_by, input$pca_color_palette, input$pca_dims)
    
    validate(need(
      input$pca_color_palette %in% names(PALETTES$all_palettes),
      "Invalid palette"
    ))
    
    tryCatch({
      p <- plot_pca(
        pca_data(),
        input$pca_color_by,
        input$pca_dims,
        PALETTES$all_palettes[[input$pca_color_palette]],
        out_dir
      )
      
      plotly::ggplotly(p, tooltip = "text")
      
    }, error = function(e) {
      shiny::validate(
        shiny::need(FALSE, paste0("Error rendering PCA plot: ", e$message))
      )
    })
  })
  
  output$pca_download_png <- downloadHandler(
    filename = function() paste0("pca_plot_", Sys.Date(), ".png"),
    content = function(file) {
      src <- file.path(out_dir, paste0("pca_plot_", input$pca_dims, "_", Sys.Date(), ".png"))
      validate(need(file.exists(src), "PNG file not ready"))
      file.copy(src, file)
    }
  )
  
  output$pca_download_pdf <- downloadHandler(
    filename = function() paste0("pca_plot_", Sys.Date(), ".pdf"),
    content = function(file) {
      src <- file.path(out_dir, paste0("pca_plot_", input$pca_dims, "_", Sys.Date(), ".pdf"))
      validate(need(file.exists(src), "PDF file not ready"))
      file.copy(src, file)
    }
  )
}