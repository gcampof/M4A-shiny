library(shinyjs)
library(dplyr)
library(archive)
suppressPackageStartupMessages(library(minfi))
# shinyOptions(progress.style="old")
# source("modules/css/load_data_css.R")
source("modules/load_data_helper.R")
source("modules/load_methylation_shiny.R")


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
        class = "p-4 border rounded shadow-sm bg-light",
        
        # Title
        h3("Quality Control - Detection P-value Analysis", class = "text-center mb-4"),
        
        # Radio buttons for detection P threshold
        div(
          class = "mb-3",
          radioButtons(
            ns("ld_qc_threshold"), 
            "Probe detection-P threshold:", 
            choices = c(
              "0.01 (strict – keeps only highest-quality probes)" = 0.01,
              "0.05 (default – widely used)" = 0.05,
              "0.10 (relaxed – allows more probes)" = 0.10
            ),
            selected = 0.05
          )
        ),
        
        # Plotly output
        fluidRow(
          column(
            width = 8,
            plotlyOutput(ns("ld_qc_barplot"), height = "400px")
          ),
          column(
            width = 4,
            div(
              class = "border rounded p-3 bg-white shadow-sm",
              h5("QC summary", class = "text-center"),
              textOutput(ns("qc_threshold_info")),
              textOutput(ns("qc_failed_summary")),
              textOutput(ns("qc_failed_samples"))
            )
          )
        ),
        br(),
        
        # Progress text
        textOutput(ns("qc_progress")),
        
        # Continue button
        div(
          class = "text-center",
          actionButton(
            ns("qc_continue"),
            "Continue",
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
load_data_server <- function(id, notify_ready) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    data_dir <- file.path(getwd(), "data")
    if (!dir.exists(data_dir)) dir.create(data_dir)
    filter_dir <- reactiveVal(file.path(getwd(), "common_files", "filter"))
    
    # Reactive values
    analysis_dir <- reactiveVal(NULL) 
    input_dir <- reactiveVal(NULL)
    results_dir <- reactiveVal(NULL)
    preprocessing_dir <- reactiveVal(NULL)
    qc_dir <- reactiveVal(NULL)
    beta_dir <- reactiveVal(NULL)
    samples_df <- reactiveVal(data.frame())
    type_selected <- reactiveVal(NULL)
    alert_message <- reactiveVal(NULL)
    qc_results <- reactiveVal(list())
    current_array_index <- reactiveVal(1)
    
    # Store the current array information that is being visualized
    array_names <- reactive({
      req(qc_results())
      names(qc_results()$rgsets)
    })
    
    # Calculate failed probes realtime, only used for visualization
    qc_failure_stats <- reactive({
      req(qc_results(), input$ld_qc_threshold)
      
      array <- array_names()[[current_array_index()]]
      detP  <- qc_results()$detections[[array]]
      thr   <- as.numeric(input$ld_qc_threshold)
      
      failed_perc <- colSums(detP > thr) / nrow(detP) * 100
      keep <- colMeans(detP) < thr
      
      list(
        threshold = thr,
        failed_perc = failed_perc,
        keep = keep,
        n_failed_samples = sum(!keep)
      )
    })
    

    
    
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
    
    
    # --- QC ADDITIONAL INFO
    output$qc_threshold_info <- renderText({
      paste("Detection P-value threshold:", qc_failure_stats()$threshold)
    })
    output$qc_failed_summary <- renderText({
      fp <- qc_failure_stats()$failed_perc
      sprintf(
        "Failed probes (%%): mean = %.2f | max = %.2f",
        mean(fp), max(fp)
      )
    })
    output$qc_failed_samples <- renderText({
      paste(
        "Samples failing QC:",
        qc_failure_stats()$n_failed_samples
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
            incProgress(0.5, detail = "Generate QC plots...")
            # generate_beta_qc(beta_and_targets$beta_path, 
            #                  beta_and_targets$targets_path, beta_dir())
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
    
    
    # --- RESET QC ARRAY INDEX ---
    observeEvent(input$run_qc, {
      current_array_index(1)
    })
    
    # --- RUN QC ---
    observeEvent(input$run_qc, {
      # update text from loading view to generating 
      shinyjs::hide("ld_idats_view")
      shinyjs::show("ld_loading_view")
      withProgress(message = "Running Quality Control", {
        incProgress(0.3, detail = "Separating unused IDATs...")
        tryCatch({
          selected_idats <- input$idat_table_rows_selected
          separate_unselected_idats(samples_df(), selected_idats, preprocessing_dir())
          parse_samplesheets(input_dir(), preprocessing_dir())
          
          incProgress(0.5, detail = "Generating Quality Control metrics for all arrays...")
          qc_results(load_qc_data_for_arrays(input$cores_selected, preprocessing_dir(), qc_dir()))
          
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
    
    
    # --- QC BARPLOT ---
    output$ld_qc_barplot <- renderPlotly({
      req(qc_results())
      req(input$ld_qc_threshold)
      
      arrays <- array_names()
      idx <- current_array_index()
      
      req(idx <= length(arrays))
      
      array <- arrays[[idx]]
      rgSet <- qc_results()$rgsets[[array]]
      detP  <- qc_results()$detections[[array]]
      
      # Build ggplot
      p <- generateDetectionPBarplot(
        array      = array,
        rgSet      = rgSet,
        detP       = detP,
        col_vector = col_vector,
        interactive = TRUE,
        threshold  = as.numeric(input$ld_qc_threshold)
      )
      
      ggplotly(p)
    })
    
    
    # --- QC PROGRESS ---
    output$qc_progress <- renderText({
      req(array_names())
      paste(
        "Array",
        current_array_index(),
        "of",
        length(array_names()),
        ":",
        array_names()[[current_array_index()]]
      )
    })
    
    # --- QC CONTINUE ---
    observeEvent(input$qc_continue, {
      idx <- current_array_index()
      arrays <- array_names()
      array <- arrays[[idx]]
      array_qc_dir <- file.path(qc_dir(), array)
      
      thr <- as.numeric(input$ld_qc_threshold)
      detP  <- qc_results()$detections[[array]]
      
      # Export failed probes per sample
      failed_perc <- colSums(detP > thr) / nrow(detP) * 100
      csv_file <- file.path(
        array_qc_dir,
        sprintf(
          "1.1-Percentage_of_failed_probes_by_sample_detection_p_%.2f_%s.csv",
          thr, array
        )
      )
      write.csv(failed_perc, file = csv_file, quote = FALSE, row.names = TRUE)
      
      # Store the threshold selected by the user
      tmp <- qc_results()
      tmp$threshold_selected[[array]] <- thr
      qc_results(tmp)
      
      # Move to the next array
      if (idx < length(arrays)) {
        current_array_index(idx + 1)
      } else {
        shinyjs::hide("ld_qc_view")
        shinyjs::show("ld_loading_view")
        
        # Generate beta matrix for all arrays
        tryCatch({
          arrays <- array_names()
          n_arrays <- length(arrays)
          processed_results <- vector("list", n_arrays)
          names(processed_results) <- arrays
          
            for (i in seq_along(arrays)) {
              array <- arrays[[i]]
              processed_results[[array]] <- generate_beta_matrix(
                cores        = input$cores_selected,
                array        = array,
                rgSet        = qc_results()$rgsets[[array]],
                detP         = qc_results()$detections[[array]],
                norm_method  = input$normalization,
                threshold    = qc_results()$threshold_selected[[array]],
                filter_dir   = filter_dir(),
                beta_dir     = beta_dir()
              )
            }
          
          # Create merged dir
          beta_merge_dir <- create_dir(file.path(beta_dir(), "merged"))
          
          # Merge SampleSheets
          merge_samplesheets(arrays, preprocessing_dir(), beta_merge_dir)

          # Merge matrices
          merge_beta_matrix(processed_results, beta_merge_dir)
          
          message("DONE!")
        }, error = function(e) {
          alert_message(list(
            type = "error",
            text = paste("Generating Beta matrix failed:", e$message)
          ))
        })
      }
    })
    
    return(list(type_selected = type_selected))
  })
}