source("modules/primary_analysis/mds/mds_utils.R")
mds_server <- function(beta_merged, targets_merged, PALETTES, out_dir, input, output, session) {
  
  mds_data <- reactive({
    req(beta_merged, targets_merged, input$mds_id_col, input$mds_top_cpgs)
    validate(need(!is.null(beta_merged), "Beta data missing"))
    validate(need(!is.null(targets_merged), "Targets data missing"))
    
    tryCatch({
      prepare_mds_data(beta_merged, targets_merged, input$mds_id_col, input$mds_top_cpgs)
    }, error = function(e) {
      shiny::validate(shiny::need(FALSE, paste0("Error preparing MDS data: ", e$message)))
      NULL
    })
  })
  
  output$mds_plot <- renderPlotly({
    req(mds_data(), input$mds_color_by, input$mds_color_palette)
    validate(need(input$mds_color_palette %in% names(PALETTES$all_palettes), "Invalid palette"))
    
    tryCatch({
      p <- plot_mds(mds_data(), input$mds_color_by, PALETTES$all_palettes[[input$mds_color_palette]], out_dir)
      plotly::ggplotly(p, tooltip = "text")
    }, error = function(e) {
      shiny::validate(shiny::need(FALSE, paste0("Error: ", e$message)))
    })
  })
  
  output$mds_download_png <- downloadHandler(
    filename = function() paste0("mds_plot_", Sys.Date(), ".png"),
    content = function(file) {
      src <- file.path(out_dir, paste0("mds_plot_", Sys.Date(), ".png"))
      validate(need(file.exists(src), "PNG file not ready"))
      file.copy(src, file)
    }
  )
  
  output$mds_download_pdf <- downloadHandler(
    filename = function() paste0("mds_plot_", Sys.Date(), ".pdf"),
    content = function(file) {
      src <- file.path(out_dir, paste0("mds_plot_", Sys.Date(), ".pdf"))
      validate(need(file.exists(src), "PDF file not ready"))
      file.copy(src, file)
    }
  )
}