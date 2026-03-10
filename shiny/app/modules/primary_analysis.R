library(shiny)
library(shinyjs)
library(bslib)
source("modules/primary_analysis/beta_matrix.R")
source("modules/primary_analysis/mds_pca.R")
source("modules/primary_analysis/umap.R")
source("modules/primary_analysis/heatmap.R")
source("modules/primary_analysis/annotations.R")




# 6. ENTRADA PARA HEATMAP
# 7. ENTRADA PARA UMAP WALKTRAP
# 8. ENTRADA PARA GLOBAL METHYLATION
# TODO: 1. PON EL PLOT DE QC DE BETA MATRIX IF TYPE_SELCTED = IDATS 
# 2. QC DEBE ESTAR A LA DERECH DE LOS BOXPLOT, COPIA ESTETICA DE LOS PLOTS POR DEFECTO



primary_analysis_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    div(
      class = "d-flex",
      style = "height: 100vh;",
      
      # --- SIDEBAR NAVIGATION ---
      div(
        id = ns("sidebar"),
        class = "bg-light border-right p-3",
        style = "width: 250px; overflow-y: auto; box-shadow: 0 0 10px rgba(0,0,0,0.1);",
        
        h4("Primary Analysis", class = "mb-4 text-center"),
        
        # Navigation buttons
        div(
          class = "d-flex flex-column gap-2",
          
          actionButton(
            ns("nav_beta_matrix"),
            "Beta Matrix",
            class = "btn btn-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_qc"),
            "QC Report",
            class = "btn btn-outline-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_mds"),
            "Multidimensional Scaling",
            class = "btn btn-outline-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_pca"),
            "Principal Component Analyisis",
            class = "btn btn-outline-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_umap"),
            "UMAP",
            class = "btn btn-outline-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_heatmap"),
            "Heatmap",
            class = "btn btn-outline-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_global"),
            "Global Methylation",
            class = "btn btn-outline-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_report"),
            "Download Report",
            class = "btn btn-outline-success w-100 text-start",
            style = "justify-content: flex-start;"
          )
        )
      ),
      
      # --- MAIN CONTENT AREA ---
      div(
        class = "flex-grow-1 p-4",
        style = "overflow-y: auto; background-color: #f8f9fa;",
        
        # Title
        div(
          style = "margin-bottom: 30px;",
          h2(textOutput(ns("view_title")), class = "mb-3"),
          hr()
        ),
        
        # --- BETA MATRIX VIEW ---
        # TODO: BUSCAR SITIO PARA MOSTRAR EL PLOT 3-pOST-FILTERED DATA QC FOR ARRAY
        # guardar como png
        shinyjs::hidden(
          div(
            id = ns("view_beta_matrix"),
            class = "content-section",
            
            uiOutput(ns("beta_matrix_tabs")),
          )
        ),
        
        # --- QC REPORT VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_qc"),
            class = "content-section",
            
            uiOutput(ns("qc_pdf_tabs")),
          )
        ),
        
        # --- MDS VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_mds"),
            class = "content-section",
            div(
              class = "card p-4 mt-3",
              fluidRow(
                column(
                  6,
                  sliderInput(
                    ns("mds_top_cpgs"),
                    "Top variable CpGs (MAD):",
                    min = 1000, max = 30000, value = 10000, step = 1000
                  )
                ),
                column(
                  6,
                  selectInput(
                    ns("mds_color_by"),
                    "Color by metadata:",
                    choices = NULL
                  )
                )
              )
            ),
            
            div(
              class = "card p-4 mt-3",
              plotlyOutput(ns("mds_plot"), height = "500px")
            )
          )
        ),
        
        
        # --- PCA VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_pca"),
            class = "content-section",
            div(
              class = "card p-4 mt-3",
              fluidRow(
                column(
                  4,
                  sliderInput(
                    ns("pca_top_cpgs"),
                    "Top variable CpGs (MAD):",
                    min = 1000, max = 30000, value = 10000, step = 1000
                  )
                ),
                column(
                  4,
                  selectInput(
                    ns("pca_color_by"),
                    "Color by metadata:",
                    choices = NULL
                  )
                ),
                column(
                  4,
                  selectInput(
                    ns("pca_dims"),
                    "PCA dimensions:",
                    choices = c("PC1 vs PC2", "PC1 vs PC3", "PC2 vs PC3"),
                    selected = "PC1 vs PC2"
                  )
                )
              )
            ),
            
            div(
              class = "card p-4 mt-3",
              plotlyOutput(ns("pca_plot"), height = "500px")
            )
          )
        ),
        
        
        # --- UMAP VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_umap"),
            class = "content-section",
            
            # ---- Controls ----
            div(
              class = "card p-4 mt-3",
              fluidRow(
                column(
                  3,
                  sliderInput(
                    ns("umap_top_cpgs"),
                    "Top variable CpGs (MAD):",
                    min = 1000, max = 30000, value = 10000, step = 1000
                  )
                ),
                column(
                  3,
                  sliderInput(
                    ns("umap_min_dist"),
                    "UMAP min_dist:",
                    min = 0.01, max = 0.99, value = 0.1, step = 0.01
                  )
                ),
                column(
                  3,
                  sliderInput(
                    ns("umap_n_neighbors"),
                    "UMAP n_neighbors:",
                    min = 2, max = 50, value = 15, step = 1
                  )
                ),
                column(
                  3,
                  selectInput(
                    ns("umap_metric"),
                    "UMAP metric:",
                    choices = c("euclidean", "manhattan", "cosine")
                  )
                ),
                column(
                  3,
                  selectInput(
                    ns("umap_id_col"),
                    "Sample ID column:",
                    choices = NULL
                  )
                ),
                column(
                  3,
                  selectInput(
                    ns("umap_color_by"),
                    "Color by metadata:",
                    choices = NULL
                  )
                ),
              ),
              fluidRow(
                column(
                  4,
                  checkboxInput(
                    ns("umap_show_labels"),
                    "Show sample labels",
                    value = FALSE
                  )
                ),
                column(
                  4,
                  checkboxInput(
                    ns("umap_show_summary"),
                    "Show parameter summary",
                    value = TRUE
                  )
                )
              )
            ),
            div(
              class = "card p-4 mt-3",
              plotlyOutput(ns("umap_plot"), height = "500px")
            ),
            
            div(
              class = "mt-3 text-center",
              downloadButton(
                ns("umap_export_pdf"),
                "Export PDF"
              ),
              downloadButton(
                ns("umap_export_png"),
                "Export PNG"
              )
            )
          )
        ),
        
        # --- HEATMAP VIEW ----
        shinyjs::hidden(
          div(
            id = ns("view_heatmap"),
            class = "content-section",
            
            tabsetPanel(
              id = ns("heatmap_main_tabs"),
              type = "tabs",
              
              # TAB 1 - HEATMAP PLOT
              tabPanel(
                title = "Heatmap Plot",
                
                # ---- Controls ----
                div(
                  class = "card p-4 mt-3",
                  
                  # Row 1 - General parameters
                  h6("General Parameters", class = "text-muted fw-bold mb-2"),
                  fluidRow(
                    column(3,
                           sliderInput(ns("heatmap_top"), "Top CpGs:",
                                       min = 100, max = 30000, value = 250, step = 100)
                    ),
                    column(3,
                           sliderInput(ns("heatmap_row_k"), "Row K:",
                                       min = 1, max = 15, value = 5, step = 1)
                    ),
                    column(3,
                           sliderInput(ns("heatmap_col_k"), "Col K:",
                                       min = 1, max = 15, value = 7, step = 1)
                    ),
                    column(3,
                           selectInput(ns("heatmap_method"), "Distance Method:",
                                       choices = c("pearson", "spearman"))
                    )
                  ),
                  
                  # Row 2 - Column selectors
                  fluidRow(
                    column(3,
                           selectInput(ns("heatmap_id_col"), "Sample ID column:",
                                       choices = NULL)
                    ),
                    column(3,
                           selectInput(ns("heatmap_sample_group_col"), "Sample group column:",
                                       choices = NULL)
                    ),
                    column(3,
                           selectInput(ns("heatmap_annotations"), "Annotations:",
                                       choices = c("Sample_Group", "Sentrix_ID", "None"),
                                       multiple = TRUE,
                                       selected = "Sample_Group")
                    ),
                    column(3,
                           div(class = "mt-4",
                               checkboxInput(ns("heatmap_show_row_names"), "Show row names", value = FALSE),
                               checkboxInput(ns("heatmap_show_col_names"), "Show col names", value = FALSE)
                           )
                    )
                  ),
                  
                  hr(),
                  
                  # Row 3 - Consensus Clustering parameters
                  h6("Consensus Clustering Parameters", class = "text-muted fw-bold mb-2"),
                  fluidRow(
                    column(3,
                           sliderInput(ns("heatmap_cc_kmax"), "Max K (cc_kmax):",
                                       min = 2, max = 20, value = 9, step = 1)
                    ),
                    column(3,
                           sliderInput(ns("heatmap_cc_reps"), "Repetitions (cc_reps):",
                                       min = 50, max = 2000, value = 500, step = 50)
                    ),
                    column(3,
                           sliderInput(ns("heatmap_cc_pItem"), "Item sampling (cc_pItem):",
                                       min = 0.5, max = 1, value = 0.8, step = 0.05)
                    ),
                    column(3,
                           sliderInput(ns("heatmap_cc_pFeature"), "Feature sampling (cc_pFeature):",
                                       min = 0.5, max = 1, value = 1, step = 0.05)
                    )
                  )
                ),
                
                # ---- Plot ----
                div(
                  class = "card p-4 mt-3",
                  plotOutput(ns("heatmap_plot"), height = "600px")
                ),
                
                # ---- Downloads ----
                div(
                  class = "mt-3 text-center",
                  downloadButton(ns("heatmap_export_pdf"), "Export HM PDF"),
                  downloadButton(ns("heatmap_export_png"), "Export HM PNG"),
                  downloadButton(ns("heatmap_export_row_clust_tsv"), "Export Row Clusters TSV")
                )
              ),
              
              # TAB 2 - CONSENSUS CLUSTER PDF viewer
              tabPanel(
                title = "Consensus Cluster",
                
                # CC-specific controls (col K selector to pick which K to display)
                div(
                  class = "card p-4 mt-3",
                  h6("Consensus Clustering Parameters", class = "text-muted fw-bold mb-2"),
                  fluidRow(
                    column(3,
                           sliderInput(ns("cc_view_col_k"), "View K:",
                                       min = 2, max = 15, value = 7, step = 1)
                    ),
                    column(3,
                           sliderInput(ns("cc_view_cc_kmax"), "Max K (cc_kmax):",
                                       min = 2, max = 20, value = 9, step = 1)
                    ),
                    column(3,
                           sliderInput(ns("cc_view_cc_reps"), "Repetitions (cc_reps):",
                                       min = 50, max = 2000, value = 500, step = 50)
                    ),
                    column(3,
                           sliderInput(ns("cc_view_cc_pItem"), "Item sampling (cc_pItem):",
                                       min = 0.5, max = 1, value = 0.8, step = 0.05)
                    )
                  )
                ),
                
                div(
                  class = "card p-4 mt-3",
                  htmlOutput(outputId = ns("consensus_cluster_pdf_tabs"), height = "600px")
                )
              )
            )
          )
        ),
        
        # --- GLOBAL METHYLATION VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_global"),
            class = "content-section",
            
            p("Global Methylation Level Analysis", class = "text-muted"),
            
            div(
              class = "row mt-4",
              column(
                width = 6,
                div(
                  class = "card p-4",
                  h5("Mean Methylation by Sample"),
                  plotOutput(ns("global_boxplot"), height = "350px")
                )
              ),
              column(
                width = 6,
                div(
                  class = "card p-4",
                  h5("Methylation Distribution"),
                  plotOutput(ns("global_density"), height = "350px")
                )
              )
            ),
            
            br(),
            actionButton(ns("action_global_export"), "Export Analysis", class = "btn btn-primary")
          )
        ),
        
        
        # --- DOWNLOAD REPORT VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_report"),
            class = "content-section",
            
            p("Generate and Download Analysis Report", class = "text-muted"),
            
            div(
              class = "card p-4 mt-3",
              h5("Report Options"),
              
              fluidRow(
                column(
                  width = 12,
                  checkboxGroupInput(
                    ns("report_sections"),
                    "Select sections to include:",
                    choices = c(
                      "QC Report" = "qc",
                      "Beta Matrix Summary" = "beta",
                      "Clustering Analysis" = "clustering",
                      "Heatmap" = "heatmap",
                      "Global Methylation" = "global",
                      "UMAP" = "umap"
                    ),
                    selected = c("qc", "beta", "global"),
                    inline = FALSE
                  )
                )
              ),
              
              br(),
              
              fluidRow(
                column(
                  width = 12,
                  textInput(
                    ns("report_title"),
                    "Report Title:",
                    value = "Methylation Analysis Report"
                  )
                )
              )
            ),
            
            br(),
            actionButton(ns("action_report_generate"), "Generate Report", class = "btn btn-success btn-lg"),
            downloadButton(ns("action_report_download"), "Download Report", class = "btn btn-success btn-lg")
          )
        )
      )
    )
  )
}

# --- SERVER ---
primary_analysis_server <- function(id, load_data_return) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # DEBUG -- HARD CODED VARIABLES
    analysis_dir <- "~/project/shiny/data/analysis_20260204_151107"
    preprocessing_dir <- "~/project/shiny/data/analysis_20260204_151107/results/1-preprocessing"
    results_dir <- "~/project/shiny/data/analysis_20260204_151107/results"
    qc_dir <- "~/project/shiny/data/analysis_20260204_151107/results/2-qc"
    beta_dir <- "~/project/shiny/data/analysis_20260204_151107/results/3-beta"
    array_names <- c("450K", "EPIC", "EPIC_V2")
    qc_results = list()
    
    tmp <- list()
    for (arr in array_names) {
      tmp[[arr]] <- readRDS(
        file.path(beta_dir, arr, paste0("001_beta_", arr, ".rds"))
      )
    }
    beta_list <- tmp
    beta_merged <- read.csv(file.path(beta_dir, "merged", "beta_merged.csv"))
    targets_merged <- read.csv(file.path(beta_dir, "merged", "targets_merged.csv"))
    
    # #CORRECT VARIABLES
    # analysis_dir <- load_data_return$analysis_dir
    # results_dir = load_data_return$results_dir
    # qc_dir = load_data_return$qc_dir
    # beta_dir = load_data_return$beta_dir
    # array_names = load_data_return$array_names
    # qc_results = load_data_return$qc_results
    # beta_list = load_data_return$beta_list
    # beta_merged = load_data_return$beta_merged
    # targets_list = load_data_return$targets_list
    # targets_merged = load_data_return$targets_merged
    
    # --- REACTIVE VALUES ---
    # Directories Primary Analysis
    beta_plots <- reactiveVal(list())
    mds_pca_dir <- reactiveVal(NULL)
    umap_dir <- reactiveVal(NULL)
    heatmap_dir <- reactiveVal(NULL)
    global_met_dir <- reactiveVal(NULL)
    
    observe({
      # Create new directories
      mds_pca_path <- create_dir(file.path(normalizePath(results_dir), "4-mds_pca"))
      umap_path <- create_dir(file.path(normalizePath(results_dir), "5-umap"))
      heatmap_path <- create_dir(file.path(normalizePath(results_dir), "6-heatmap"))
      global_met_path <- create_dir(file.path(normalizePath(results_dir), "7-global_methylation"))
      
      # Assign to the reactive values
      mds_pca_dir(mds_pca_path)
      umap_dir(umap_path)
      heatmap_dir(heatmap_path)
      global_met_dir(global_met_path)
      
      # Add QC dir path for pdf viewer
      addResourcePath(
        prefix = "qc_reports",
        directoryPath = normalizePath(qc_dir)
      )
      
      # Initialize metadata color choices dynamically
      meta_cols <- colnames(targets_merged)
      updateSelectInput(session, "pca_color_by", choices = meta_cols, selected = meta_cols[1])
      updateSelectInput(session, "mds_color_by", choices = meta_cols, selected = meta_cols[1])
      updateSelectInput(session, "umap_color_by", choices = meta_cols, selected = meta_cols[1])
      
      # Initializa sample ID choices dinamically (default "ID")
      updateSelectInput(session, "umap_id_col", choices = meta_cols, 
                        selected = if ("ID" %in% meta_cols) "ID" else meta_cols[1])    
      })
    
    # --- NAVIGATION BUTTONS ---
    observeEvent(input$nav_beta_matrix, {
      current_view("beta_matrix")
      update_active_button("nav_beta_matrix")
      show_view("view_beta_matrix", "Beta Matrix")
    })
    
    observeEvent(input$nav_qc, {
      current_view("qc")
      update_active_button("nav_qc")
      show_view("view_qc", "QC Report")
    })
    
    observeEvent(input$nav_mds, {
      current_view("mds")
      update_active_button("nav_mds")
      show_view("view_mds", "Multidimensional Scaling (MDS)")
    })
    
    observeEvent(input$nav_pca, {
      current_view("pca")
      update_active_button("nav_pca")
      show_view("view_pca", "Principal Component Analysis (PCA)")
    })
    
    observeEvent(input$nav_umap, {
      current_view("umap")
      update_active_button("nav_umap")
      show_view("view_umap", "UMAP")
    })
    
    observeEvent(input$nav_heatmap, {
      current_view("heatmap")
      update_active_button("nav_heatmap")
      show_view("view_heatmap", "Heatmap")
    })
    
    observeEvent(input$nav_global, {
      current_view("global")
      update_active_button("nav_global")
      show_view("view_global", "Global Methylation")
    })
    

    # --- BETA MATRIX BOXPLOT UI ---
    output$beta_matrix_tabs <- renderUI({
      tabs <- lapply(array_names, function(arr) {
        tabPanel(
          title = arr,
          shinyjs::hidden(
            div(
              id = ns(paste0("spinner_", arr)),
              class = "spinner-border text-primary",
              role = "status",
              style = "width: 5rem; height: 5rem;"
            )
          ),
          
          br(),
          plotlyOutput(
            outputId = ns(paste0("pa_beta_boxplot_", arr)),
            height = "600px"
          ),
          )
      })
      
      do.call(
        tabsetPanel,
        c(
          list(
            id = ns("beta_matrix_tabset"),
            type = "tabs"
          ),
          tabs
        )
      )
    })
    
    # --- BETA MATRIX BOXPLOT LOGIC ---
    observeEvent(input$beta_matrix_tabset, {
      arr <- input$beta_matrix_tabset
      req(arr)
      beta <- beta_list[[arr]]
      
      spinner_id <- paste0("spinner_", arr)
      shinyjs::show(spinner_id)
      
      # Generate plot if not already cached
      if (is.null(beta_plots[[arr]])) {
        beta_plots[[arr]] <- generate_beta_boxplots(arr, beta, beta_dir)
      } 
      output_id <- paste0("pa_beta_boxplot_", arr)
      output[[output_id]] <- renderPlotly({
        beta_plots[[arr]]
      })
      
      shinyjs::hide(spinner_id)
    })
    
    
    # --- QC PDF VIEWER UI ---
    output$qc_pdf_tabs <- renderUI({
      tabs <- lapply(array_names, function(arr) {
        tabPanel(
          title = arr,
          htmlOutput(outputId = ns(paste0("pa_qc_viewer", arr)),
                     height = "600px")
        )
      })
      
      do.call(
        tabsetPanel,
        c(
          list(
            id = ns("qc_pdf_tabset"),
            type = "tabs"
          ),
          tabs
        )
      )
    })
    
    # --- QC PDF VIEWER LOGIC ---
    observeEvent(input$qc_pdf_tabset, {
      arr <- input$qc_pdf_tabset
      req(arr)
      beta <- beta_list[[arr]]
      
      src <- file.path("qc_reports", arr, paste0("2.0-QC_Report",arr,".pdf"))
      
      output_id <- paste0("pa_qc_viewer", arr)
      output[[output_id]] <- renderText({
        return(paste('<iframe style="height:600px; width:100%" src="', src, '"></iframe>', sep = ""))
      })
    })
    
    # --- MDS PLOT LOGIC ---
    mds_data <- reactive({
      req(beta_merged, input$mds_top_cpgs)
      prepare_mds_data(beta_merged, input$mds_top_cpgs)
    })
    
    output$mds_plot <- renderPlotly({
      req(mds_data(), input$mds_color_by)
      plot_mds(mds_data(), targets_merged, input$mds_color_by)
    })
    
    
    # --- PCA PLOT LOGIC ---
    # Generate pca df
    pca_data <- reactive({
      req(beta_merged, input$pca_top_cpgs)
      prepare_pca_data(beta_merged, input$pca_top_cpgs)
    })
    
    # Plot pca
    output$pca_plot <- renderPlotly({
      req(pca_data(), input$pca_color_by)
      plot_pca(pca_data(), targets_merged, input$pca_color_by, input$pca_dims)
    })
    
    
    # --- UMAP PLOT LOGIC ---
    # Genereate umap df
    umap_data <- reactive({
      req(beta_merged, targets_merged, input$umap_top_cpgs, input$umap_min_dist,
          input$umap_n_neighbors, input$umap_metric, input$umap_id_col)
      
      prepare_umap_data(beta_merged, targets_merged, input$umap_top_cpgs,
                        input$umap_min_dist, input$umap_n_neighbors, input$umap_metric,
                        id_col = input$umap_id_col)
    })
    
    # Plot umap
    output$umap_plot <- renderPlotly({
      req(umap_data(), input$umap_color_by)
      plot_umap(
        umap_df = umap_data(),
        color_by = input$umap_color_by,
        show_labels = input$umap_show_labels,
        show_summary = input$umap_show_summary,
        top_cpgs = input$umap_top_cpgs,
        min_dist = input$umap_min_dist,
        n_neighbors = input$umap_n_neighbors,
        metric = input$umap_metric
      )
    })
    
    
    # --- HEATMAP PLOT LOGIC ---
    # Generate heatmap data
    heatmap_data <- reactive({
      req(beta_merged, targets_merged, input$umap_top_cpgs, input$umap_min_dist,
          input$umap_n_neighbors, input$umap_metric, input$umap_id_col)
      prepare_umap_data(beta_merged, targets_merged, input$umap_top_cpgs,
                        input$umap_min_dist, input$umap_n_neighbors, input$umap_metric,
                        id_col = input$umap_id_col)
    })
    
    # Plot heatmap
    output$heatmap_plot <- renderPlotly({
      req(umap_data(), input$umap_color_by)
      
      plot_umap(
        umap_df = umap_data(),
        color_by = input$umap_color_by,
        show_labels = input$umap_show_labels,
        show_summary = input$umap_show_summary,
        top_cpgs = input$umap_top_cpgs,
        min_dist = input$umap_min_dist,
        n_neighbors = input$umap_n_neighbors,
        metric = input$umap_metric
      )
    })
    
    
    
    # --- HELPER FUNCTIONS ---
    update_active_button <- function(active_id) {
      shinyjs::removeClass(selector = ".content-section", class = "active")
      
      all_buttons <- c("nav_beta_matrix", "nav_qc", "nav_mds", "nav_pca", 
                       "nav_clustering", "nav_heatmap", "nav_global", 
                       "nav_umap", "nav_report")
      
      for (btn in all_buttons) {
        if (btn == active_id) {
          shinyjs::addClass(id = btn, class = "btn-primary")
          shinyjs::removeClass(id = btn, class = "btn-outline-primary")
        } else {
          shinyjs::removeClass(id = btn, class = "btn-primary")
          shinyjs::addClass(id = btn, class = "btn-outline-primary")
        }
      }
    }
    
    show_view <- function(view_id, title_text) {
      all_views <- c("view_beta_matrix", "view_qc", "view_mds", "view_pca", 
                     "view_clustering", "view_heatmap", "view_global",
                     "view_umap", "view_report")
      
      for (view in all_views) {
        shinyjs::hide(view)
      }
      
      shinyjs::show(view_id)
      output$view_title <- renderText(title_text)
    }
    
    
    # --- PLACEHOLDER ACTIONS ---
    # You can replace these with actual analysis functions
    
    observeEvent(input$action_clustering_generate, {
      showNotification("Generating clustering analysis...", type = "message")
    })
    
    
    observeEvent(input$action_report_generate, {
      showNotification("Generating report...", type = "message")
    })
    
    output$action_report_download <- downloadHandler(
      filename = function() {
        paste0("methylation_report_", Sys.Date(), ".html")
      },
      content = function(file) {
        writeLines("Placeholder report content", file)
      }
    )
    
    # Initialize first view
    show_view("view_beta_matrix", "Beta Matrix")
    update_active_button("nav_beta_matrix")
    
  })
}