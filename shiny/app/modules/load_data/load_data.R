library(shinyjs)
library(dplyr)
library(archive)
library(doParallel)
library(foreach)
suppressPackageStartupMessages(library(minfi))
# shinyOptions(progress.style="old")
# source("modules/css/load_data_css.R")
source("modules/load_data/load_data_helper.R")


load_data_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    shinyjs::useShinyjs(),
    
    # --- TITLE ---
    h2("Load Your Data", class = "text-center mb-4"),
    
    # --- MESSAGE ALERTS ---
    uiOutput(ns("alert_message")),
    
    # --- Main view  ---
    div(
      style = "max-width: 900px; width: 100%;",
      id = ns("ld_main_view"),
      class = "p-4 border rounded shadow-sm bg-light",
      p("Select the type of files you want to load:", class = "text-center text-muted mb-4"),
      
      # Selection buttons centered
      div(
        class = "d-flex justify-content-center gap-4 mb-5",
        style = "flex-wrap: wrap;",
        actionButton(
          ns("load_beta"), 
          "Beta Matrix\n+ Targets file", 
          class = "btn btn-primary btn-lg",
          style = "width: 280px; height: 280px; font-size: 18px; border-radius: 8px;"
        ),
        actionButton(
          ns("load_idats"), 
          "IDATS\n+ Samplesheet",
          class = "btn btn-primary btn-lg",
          style = "width: 280px; height: 280px; font-size: 18px; border-radius: 8px;"
        )
      ),
      
      # File upload centered
      div(
        class = "d-flex justify-content-center mb-4",
        uiOutput(ns('zipfile_ui'))
      ),
      
      # Load Data button - bigger and centered with spinner
      div(
        class = "text-center",
        style = "position: relative;",
        actionButton(
          inputId = ns("confirm_load"),
          label = "Load Data",
          class = "btn btn-success btn-lg",
          style = "width: 280px; height: 60px; font-size: 18px; border-radius: 8px; position: relative;",
          disabled = TRUE
        ),
        # Spinner overlay (hidden by default)
        div(
          id = ns("load_spinner"),
          class = "spinner-border spinner-border-sm text-primary",
          role = "status",
          style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); display: none; z-index: 10;",
          span(class = "visually-hidden", "Loading...")
        )
      )
    ),
    
    # --- IDATS View ---
    shinyjs::hidden(
      div(
        id = ns("ld_idats_view"),
        class = "p-4 border rounded shadow-sm bg-light",
        h3("Detected IDAT Samples", class = "text-center mb-4"),
        
        # Select all checkbox
        div(
          class = "text-center mb-3",
          checkboxInput(ns("select_all"), "Select all samples", value = FALSE)
        ),
        
        # DataTable
        DTOutput(ns("idat_table")),
        br(),
        
        # Normalization + detection threshold + CPU selectors
        fluidRow(
          column(
            width = 4,
            selectInput(
              ns("normalization"), 
              "Normalization method:", 
              choices = c("ssnoob", "raw", "illumina", "quantile", "funnorm"), 
              selected = "ssnoob"
            )
          ),
          column(
            width = 4,
            uiOutput(ns("cpu_selector")),
            textOutput(ns("available_cores_text"))
          )
        ),
        br(),
        
        # Run QC button
        div(
          class = "text-center",
          actionButton(
            ns("run_qc"), 
            "Run QC", 
            class = "btn btn-primary btn-lg",
            style = "width: 280px; height: 60px; font-size: 18px; border-radius: 8px;"
          )
        )
      )
    ),
    
    # --- QC VIEW ---
    shinyjs::hidden(
      div(
        id = ns("ld_qc_view"),
        class = "content-section",
        uiOutput(ns("qc_threshold_tabs")),
        br(),
        div(
          class = "text-center",
          actionButton(
            ns("qc_continue"),
            "Generate Beta Matrix",
            class = "btn btn-primary btn-lg",
            style = "width: 280px; height: 60px; font-size: 18px; border-radius: 8px;"
          )
        )
      )
    ),
    
    
    # --- LOADING View ---
    shinyjs::hidden(
      div(
        id = ns("ld_loading_view"),
        class = "text-center p-5",
        
        # Bootstrap spinner
        div(
          class = "spinner-border text-primary",
          role = "status",
          style = "width: 5rem; height: 5rem;"
        ),
        # TODO: A CONSOLE THAT SHOWS THE MESSAGES THAT WE ARE PRINTING, 
        # ALSO IT SHOULD SAVE THE LOG INTO A LOG.TXT FILE
        div(
          style = "border: 1px solid #ccc; padding: 10px; height: 300px; overflow-y: scroll; background-color: #f7f7f7;",
          id = "console_log"
        ),
        br(), br(),
      )
    )
  )
}


# --- SERVER ---
load_data_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    data_dir <- file.path(getwd(), "data")
    if (!dir.exists(data_dir)) dir.create(data_dir)
    filter_dir <- reactiveVal(file.path(getwd(), "common_files", "filter"))
    
    # --- REACTIVE VALUES ---
    # Directories
    analysis_dir <- reactiveVal(NULL) 
    input_dir <- reactiveVal(NULL)
    results_dir <- reactiveVal(NULL)
    preprocessing_dir <- reactiveVal(NULL)
    qc_dir <- reactiveVal(NULL)
    beta_dir <- reactiveVal(NULL)
    samples_df <- reactiveVal(data.frame())
    type_selected <- reactiveVal(NULL)
    alert_message <- reactiveVal(NULL)
    
    # QC reactive values
    qc_results <- reactiveVal(list())
    array_names <- reactiveVal(list())

    # Beta reactive values
    beta_list <- reactiveVal(list())
    beta_merged <- reactiveVal(NULL)
    targets_list <- reactiveVal(list())
    targets_merged <- reactiveVal(NULL)
    
    # --- RENDER ZIP FILE ---
    output$zipfile_ui <- renderUI({
      fileInput(
        inputId = ns("zipfile"),
        label = "Upload your ZIP files",
        accept = ".zip",
        multiple = TRUE,
        width = "280px"
      )
    })
    
    # --- ALERT MESSAGE ---
    output$alert_message <- renderUI({
      msg <- alert_message()
      if (is.null(msg)) return(NULL)
      alert_class <- if (msg$type == "error") "alert-danger" else "alert-warning"
      div(
        class = paste("alert", alert_class, "alert-dismissible fade show"),
        role = "alert",
        msg$text,
      )
    })
    
    # --- CPU SELECTOR ---
    available_cores <- parallel::detectCores()
    output$cpu_selector <- renderUI({
      sliderInput(
        ns("cores_selected"),
        "Select number of CPUs to use:",
        min = 1,
        max = available_cores,
        value = min(available_cores, 4),
        step = 1
      )
    })
    
    output$available_cores_text <- renderText({
      paste("Detected", available_cores, "available CPU cores")
    })
    
    # --- IDAT TABLE ---
    output$idat_table <- DT::renderDataTable({
      req(samples_df())
      datatable(
        samples_df(),
        selection = list(mode = "multiple"),
        filter = "top",
        options = list(pageLength = 10)
      )
    })
    
    # --- TYPE SELECTOR ---
    observeEvent(input$load_beta, {
      type_selected("BETA")
      alert_message(NULL)
      shinyjs::removeClass("load_beta", "btn-primary")
      shinyjs::addClass("load_beta", "btn-success")
      shinyjs::removeClass("load_idats", "btn-success")
      shinyjs::addClass("load_idats", "btn-primary")
    })
    
    observeEvent(input$load_idats, {
      type_selected("IDATS")
      alert_message(NULL)
      shinyjs::removeClass("load_idats", "btn-primary")
      shinyjs::addClass("load_idats", "btn-success")
      shinyjs::removeClass("load_beta", "btn-success")
      shinyjs::addClass("load_beta", "btn-primary")
    })
    
    # --- ENABLE LOAD DATA BUTTON ---
    observe({
      type_ok <- !is.null(type_selected())
      zip_ok <- !is.null(input$zipfile) && nrow(input$zipfile) > 0
      
      if (type_ok && zip_ok) {
        shinyjs::enable("confirm_load")
        shinyjs::addClass("confirm_load", "btn-success")
        shinyjs::removeClass("confirm_load", "btn-secondary")
      } else {
        shinyjs::disable("confirm_load")
        shinyjs::removeClass("confirm_load", "btn-success")
        shinyjs::addClass("confirm_load", "btn-secondary")
      }
    })
    
    # --- FILE UPLOAD ---
    observeEvent(input$confirm_load, {
      # Show loading view
      shinyjs::hide("ld_main_view")
      shinyjs::show("ld_loading_view")
      
      withProgress(message = "Loading data", {
        
        incProgress(0.1, detail = "Preparing directories...")
        
        # Create all required folders
        current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
        
        incProgress(0.2, detail = "Creating directory structure...")
        filter_path <- file.path()
        analysis_path <- file.path(data_dir, paste0("analysis_", current_time))
        analysis_path <- create_dir(analysis_path)
        input_path <- create_dir(file.path(analysis_path, "input"))
        results_path <- create_dir(file.path(analysis_path, "results"))
        preprocessing_path <- create_dir(file.path(results_path, "1-preprocessing"))
        qc_path <- create_dir(file.path(results_path, "2-qc"))
        beta_path <- create_dir(file.path(results_path, "3-beta"))
        
        # Assign to the reactive values
        analysis_dir(analysis_path)
        input_dir(input_path)
        results_dir(results_path)
        preprocessing_dir(preprocessing_path)
        qc_dir(qc_path)
        beta_dir(beta_path)
        
        tryCatch({
          incProgress(0.5, detail = "Extracting files from zip...")
          zip_paths <- input$zipfile$datapath
          archive::archive_extract(zip_paths, dir = input_dir())
          
          if (type_selected() == "IDATS") {
            incProgress(0.9, detail = "Generating Sample IDAT Summary..")
            idats_dir <- parse_idat_files(input_dir(), preprocessing_dir())
            order_idat_per_array(idats_dir, preprocessing_dir())
            temp_df <- generate_idat_dataframe(preprocessing_dir())
            samples_df(temp_df)
            
            # Show IDATS view
            incProgress(1.0, detail = "Complete!")
          } else if (type_selected() == "BETA") {
            incProgress(0.5, detail = "Extracting BETA files...")
            beta_and_targets <- extract_beta_and_targets(input_dir(), beta_dir())
            beta_list(beta_and_targets$beta)
            targets_list(beta_and_targets$targets)
            incProgress(1.0, detail = "Complete!")
          }
        }, error = function(e) {
          # Remove created analysis directory if it exists
          # TODO: Leave it commented for testing
          # if (dir.exists(analysis_dir)) {
          #   unlink(analysis_dir, recursive = TRUE, force = TRUE)
          # }
          
          # Reset the ZIP upload input 
          output$zipfile_ui <- renderUI({
            fileInput(
              inputId = ns("zipfile"),
              label = "Upload your ZIP files",
              accept = ".zip",
              multiple = TRUE,
              width = "280px"
            )
          })          
          alert_message(list(
            type = "error",
            text = paste("Processing failed:", e$message)
          ))
        })
      })
      
      # Change views
      if (type_selected() == "IDATS") {
        shinyjs::hide("ld_loading_view")
        shinyjs::show("ld_idats_view")
      }
    })
    
    # --- SELECT ALL ---
    observeEvent(input$select_all, {
      df <- samples_df()
      if (nrow(df) == 0) return()
      
      proxy <- DT::dataTableProxy("idat_table", session)
      
      if (isTRUE(input$select_all)) {
        proxy %>% DT::selectRows(seq_len(nrow(df)))
      } else {
        proxy %>% DT::selectRows(NULL)
      }
    })
    
    
    # --- RUN QC ---
    observeEvent(input$run_qc, {
      shinyjs::hide("ld_idats_view")
      shinyjs::show("ld_loading_view")
      withProgress(message = "Running Quality Control", {
        incProgress(0.3, detail = "Separating unused IDATs...")
        
        tryCatch({
          selected_idats <- input$idat_table_rows_selected
          separate_unselected_idats(samples_df(), selected_idats, preprocessing_dir())
          parse_samplesheets(input_dir(), preprocessing_dir())
          
          incProgress(0.5, detail = "Generating Quality Control metrics for all arrays...")
          
          qc_res_temp <- load_qc_data_for_arrays(input$cores_selected, preprocessing_dir(), qc_dir())
          qc_results(qc_res_temp$qc_results)
          array_names(qc_res_temp$arrays_used)
          
          incProgress(1.0, detail = "Complete!")
          shinyjs::hide("ld_loading_view")
          shinyjs::show("ld_qc_view")
        }, error = function(e) {
          alert_message(list(
            type = "error",
            text = paste("Running QC failed:", e$message)
          ))
        })
      })
    })
    
    
    
    
    # --- QC TABS VIEW UI ---
    output$qc_threshold_tabs <- renderUI({
      req(array_names())
      
      tabs <- lapply(array_names(), function(arr) {
        tabPanel(
          title = arr,
          br(),
          
          radioButtons(
            ns(paste0("ld_qc_threshold_", arr)),
            "Probe detection-P threshold:",
            choices = c(
              "0.01 (strict)" = 0.01,
              "0.05 (default)" = 0.05,
              "0.10 (relaxed)" = 0.10
            ),
            selected = 0.05,
            inline = TRUE
          ),
          
          fluidRow(
            column(
              width = 9,
              plotlyOutput(
                outputId = ns(paste0("ld_qc_barplot_", arr)),
                height = "600px"
              )
            ),
            column(
              width = 3,
              uiOutput(ns(paste0("ld_qc_stats_", arr)))
            )
          )
        )
      })
      
      do.call(
        tabsetPanel,
        c(
          list(
            id = ns("qc_threshold_tabset"),
            type = "tabs"
          ),
          tabs
        )
      )
    })
    
    # -- QC TABS VIEW LOGIC ---
    observe({
      req(qc_results(), array_names())
      
      lapply(array_names(), function(arr) {
        
        local({
          array    <- arr
          plot_id  <- paste0("ld_qc_barplot_", array)
          stats_id <- paste0("ld_qc_stats_", array)
          thr_id   <- paste0("ld_qc_threshold_", array)
          
          # ---- Plot ----
          output[[plot_id]] <- renderPlotly({
            req(input[[thr_id]])
            
            detP  <- qc_results()$detections[[array]]
            rgSet <- qc_results()$rgsets[[array]]
            thr   <- as.numeric(input[[thr_id]])
            
            p <- generate_detection_p_barplot(
              array     = array,
              rgSet     = rgSet,
              detP      = detP,
              threshold = thr
            )
            
            ggplotly(p)
          })
          
          # ---- QC stats box ----
          output[[stats_id]] <- renderUI({
            req(input[[thr_id]])
            
            detP <- qc_results()$detections[[array]]
            thr  <- as.numeric(input[[thr_id]])
            
            failed_perc <- colSums(detP > thr) / nrow(detP) * 100
            keep        <- colMeans(detP) < thr
            
            tagList(
              div(
                class = "border rounded p-3 bg-light",
                
                h5("QC summary"),
                
                p(strong("Threshold: "), thr),
                p(
                  strong("Mean failed probes (%): "),
                  sprintf("%.2f", mean(failed_perc))
                ),
                p(
                  strong("Max failed probes (%): "),
                  sprintf("%.2f", max(failed_perc))
                ),
                p(
                  strong("Samples failing QC: "),
                  sum(!keep)
                )
              )
            )
          })
        })
      })
    })
    
    
    # --- QC CONTINUE ---
    observeEvent(input$qc_continue, {
      # ---- Move forward ----
      shinyjs::hide("ld_qc_view")
      shinyjs::show("ld_loading_view")
      
      arrays <- array_names()
      qc_res <- qc_results()
      
      # Save output
      n_arrays <- length(arrays)
      processed_results <- vector("list", n_arrays)
      names(processed_results) <- arrays      
      
      tryCatch({
        
        for (array in arrays) {
          message("Saving QC results for array: ", array)
          array_qc_dir <- file.path(qc_dir(), array)
          
          # ---- Get threshold selected for THIS array ----
          thr_id <- paste0("ld_qc_threshold_", array)
          thr <- as.numeric(input[[thr_id]])
          
          detP  <- qc_res$detections[[array]]
          rgSet <- qc_res$rgsets[[array]]
          
          # Save Percentage of failed probes
          failed_perc <- colSums(detP > thr) / nrow(detP) * 100
          
          csv_file <- file.path(
            array_qc_dir,
            sprintf(
              "1.1-Percentage_of_failed_probes_by_sample_detection_p_%.2f_%s.csv",
              thr, array
            )
          )
          write.csv(failed_perc, file = csv_file, quote = FALSE, row.names = TRUE)
          
          # ---- Save QC plot ----
          p <- generate_detection_p_barplot(
            array     = array,
            rgSet     = rgSet,
            detP      = detP,
            threshold = thr
          )
          
          plot_file <- file.path(
            array_qc_dir,
            sprintf(
              "1.0-Detection_P_barplot_threshold_%.2f_%s.png",
              thr, array
            )
          )
          ggsave(filename = plot_file, plot = p, width = 12, height = 6, dpi = 300)
          
          # ---- Store threshold ----
          qc_res$threshold_selected[[array]] <- thr
          
          
          # --- Generate Beta matrix ---
          processed_results[[array]] <- generate_beta_matrix(
            cores       = input$cores_selected,
            array       = array,
            rgSet       = qc_res$rgsets[[array]],
            detP        = qc_res$detections[[array]],
            norm_method = input$normalization,
            threshold   = qc_res$threshold_selected[[array]],
            filter_dir  = filter_dir(),
            beta_dir    = beta_dir()
          )
          
          # Store beta on list
          bl <- beta_list()
          bl[[array]] <- processed_results$beta
          beta_list(bl)
        }
        
        qc_results(qc_res)
        
        # Merge Beta matrices
        beta_merge_dir <- create_dir(file.path(beta_dir(), "merged"))
        beta_mg <- merge_beta_matrix(processed_results, beta_merge_dir)
        beta_merged(beta_mg)
        
        # Merge targets file
        targets_result <- merge_samplesheets(arrays, preprocessing_dir(), beta_merge_dir)
        targets_list(targets_result$targets_individual)
        targets_merged(targets_result$targets_merged)    
        
        message("QC and beta generation completed successfully")
        
      }, error = function(e) {
        alert_message(list(
          type = "error",
          text = paste("Generating Beta matrix failed:", e$message)
        ))
      })
    })
    
    return(list(
      analysis_dir = analysis_dir,
      results_dir = results_dir,
      qc_dir = qc_dir,
      beta_dir = beta_dir,
      qc_results = qc_results,
      array_names = array_names,
      beta_list = beta_list,
      beta_merged = beta_merged,
      targets_list = targets_list,
      targets_merged = targets_merged
    ))
  })
}