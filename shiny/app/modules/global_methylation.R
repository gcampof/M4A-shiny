library(shiny)
library(shinyjs)
library(bslib)

global_methylation_ui <- function(id) {
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
        
        h4("Global Methylation", class = "mb-4 text-center"),
        
        # Navigation buttons
        div(
          class = "d-flex flex-column gap-2",
          
          actionButton(
            ns("nav_beta_matrix"),
            "Beta Matrix",
            class = "btn btn-outline-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_qc"),
            "QC Report",
            class = "btn btn-outline-primary w-100 text-start",
            style = "justify-content: flex-start;"
          ),
          actionButton(
            ns("nav_clustering"),
            "Hierarchical Clustering",
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
            ns("nav_umap"),
            "UMAP",
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
        shinyjs::hidden(
          div(
            id = ns("view_beta_matrix"),
            class = "content-section",
            
            p("Beta Matrix visualization and analysis", class = "text-muted"),
            
            div(
              class = "row mt-4",
              column(
                width = 12,
                div(
                  class = "card p-4",
                  h5("Beta Matrix Summary"),
                  p("Dimensions: Placeholder"),
                  p("Range: 0 to 1"),
                  p("Missing values: Placeholder")
                )
              )
            ),
            
            br(),
            actionButton(ns("action_beta_export"), "Export Beta Matrix", class = "btn btn-primary")
          )
        ),
        
        # --- QC REPORT VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_qc"),
            class = "content-section",
            
            p("Quality Control Report", class = "text-muted"),
            
            div(
              class = "row mt-4",
              column(
                width = 6,
                div(
                  class = "card p-4",
                  h5("Detection P-values"),
                  plotOutput(ns("qc_detection_plot"), height = "300px")
                )
              ),
              column(
                width = 6,
                div(
                  class = "card p-4",
                  h5("Bisulfite Conversion"),
                  plotOutput(ns("qc_bisulfite_plot"), height = "300px")
                )
              )
            ),
            
            br(),
            actionButton(ns("action_qc_export"), "Export QC Report", class = "btn btn-primary")
          )
        ),
        
        # --- HIERARCHICAL CLUSTERING VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_clustering"),
            class = "content-section",
            
            p("Hierarchical Clustering Analysis", class = "text-muted"),
            
            div(
              class = "card p-4 mt-3",
              fluidRow(
                column(
                  width = 6,
                  selectInput(
                    ns("cluster_method"),
                    "Clustering Method:",
                    choices = c("Ward", "Complete", "Average", "Single")
                  )
                ),
                column(
                  width = 6,
                  sliderInput(
                    ns("cluster_samples"),
                    "Number of Top Variable Probes:",
                    min = 100,
                    max = 10000,
                    value = 5000,
                    step = 500
                  )
                )
              )
            ),
            
            div(
              class = "card p-4 mt-3",
              h5("Dendrogram"),
              plotOutput(ns("clustering_plot"), height = "500px")
            ),
            
            br(),
            actionButton(ns("action_clustering_generate"), "Generate Clustering", class = "btn btn-primary"),
            actionButton(ns("action_clustering_export"), "Export Results", class = "btn btn-secondary")
          )
        ),
        
        # --- HEATMAP VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_heatmap"),
            class = "content-section",
            
            p("Heatmap Visualization", class = "text-muted"),
            
            div(
              class = "card p-4 mt-3",
              fluidRow(
                column(
                  width = 6,
                  sliderInput(
                    ns("heatmap_probes"),
                    "Number of Probes:",
                    min = 100,
                    max = 5000,
                    value = 1000,
                    step = 100
                  )
                ),
                column(
                  width = 6,
                  radioButtons(
                    ns("heatmap_scale"),
                    "Scale method:",
                    choices = c("Row-wise", "Column-wise", "None"),
                    inline = TRUE
                  )
                )
              )
            ),
            
            div(
              class = "card p-4 mt-3",
              h5("Heatmap Preview"),
              plotOutput(ns("heatmap_plot"), height = "600px")
            ),
            
            br(),
            actionButton(ns("action_heatmap_generate"), "Generate Heatmap", class = "btn btn-primary"),
            actionButton(ns("action_heatmap_export"), "Export Heatmap", class = "btn btn-secondary")
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
        
        # --- UMAP VIEW ---
        shinyjs::hidden(
          div(
            id = ns("view_umap"),
            class = "content-section",
            
            p("UMAP Dimensionality Reduction", class = "text-muted"),
            
            div(
              class = "card p-4 mt-3",
              fluidRow(
                column(
                  width = 4,
                  sliderInput(
                    ns("umap_neighbors"),
                    "Number of Neighbors:",
                    min = 5,
                    max = 50,
                    value = 15
                  )
                ),
                column(
                  width = 4,
                  sliderInput(
                    ns("umap_min_dist"),
                    "Minimum Distance:",
                    min = 0.01,
                    max = 1,
                    value = 0.1,
                    step = 0.05
                  )
                ),
                column(
                  width = 4,
                  sliderInput(
                    ns("umap_probes"),
                    "Top Probes:",
                    min = 100,
                    max = 10000,
                    value = 5000,
                    step = 500
                  )
                )
              )
            ),
            
            div(
              class = "card p-4 mt-3",
              h5("UMAP Plot"),
              plotOutput(ns("umap_plot"), height = "500px")
            ),
            
            br(),
            actionButton(ns("action_umap_generate"), "Generate UMAP", class = "btn btn-primary"),
            actionButton(ns("action_umap_export"), "Export UMAP", class = "btn btn-secondary")
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
global_methylation_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive value to track current view
    current_view <- reactiveVal("beta_matrix")
    
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
    
    observeEvent(input$nav_clustering, {
      current_view("clustering")
      update_active_button("nav_clustering")
      show_view("view_clustering", "Hierarchical Clustering")
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
    
    observeEvent(input$nav_umap, {
      current_view("umap")
      update_active_button("nav_umap")
      show_view("view_umap", "UMAP")
    })
    
    observeEvent(input$nav_report, {
      current_view("report")
      update_active_button("nav_report")
      show_view("view_report", "Download Report")
    })
    
    # --- HELPER FUNCTIONS ---
    update_active_button <- function(active_id) {
      shinyjs::removeClass(selector = ".content-section", class = "active")
      
      all_buttons <- c("nav_beta_matrix", "nav_qc", "nav_clustering", "nav_heatmap", 
                       "nav_global", "nav_umap", "nav_report")
      
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
      all_views <- c("view_beta_matrix", "view_qc", "view_clustering", "view_heatmap", 
                     "view_global", "view_umap", "view_report")
      
      for (view in all_views) {
        shinyjs::hide(view)
      }
      
      shinyjs::show(view_id)
      output$view_title <- renderText(title_text)
    }
    
    # --- PLACEHOLDER ACTIONS ---
    # You can replace these with actual analysis functions
    
    observeEvent(input$action_beta_export, {
      showNotification("Beta Matrix exported!", type = "message")
    })
    
    observeEvent(input$action_qc_export, {
      showNotification("QC Report exported!", type = "message")
    })
    
    observeEvent(input$action_clustering_generate, {
      showNotification("Generating clustering analysis...", type = "message")
    })
    
    observeEvent(input$action_clustering_export, {
      showNotification("Clustering results exported!", type = "message")
    })
    
    observeEvent(input$action_heatmap_generate, {
      showNotification("Generating heatmap...", type = "message")
    })
    
    observeEvent(input$action_heatmap_export, {
      showNotification("Heatmap exported!", type = "message")
    })
    
    observeEvent(input$action_global_export, {
      showNotification("Global methylation analysis exported!", type = "message")
    })
    
    observeEvent(input$action_umap_generate, {
      showNotification("Generating UMAP...", type = "message")
    })
    
    observeEvent(input$action_umap_export, {
      showNotification("UMAP exported!", type = "message")
    })
    
    observeEvent(input$action_report_generate, {
      showNotification("Generating report...", type = "message")
    })
    
    output$action_report_download <- downloadHandler(
      filename = function() {
        paste0("methylation_report_", Sys.Date(), ".html")
      },
      content = function(file) {
        # TODO: Implement actual report generation
        writeLines("Placeholder report content", file)
      }
    )
    
    # --- PLACEHOLDER PLOTS ---
    output$qc_detection_plot <- renderPlot({
      plot(1:10, rnorm(10), main = "Detection P-values", xlab = "Sample", ylab = "P-value")
    })
    
    output$qc_bisulfite_plot <- renderPlot({
      plot(1:10, rnorm(10) + 1, main = "Bisulfite Conversion", xlab = "Sample", ylab = "Conversion %")
    })
    
    output$clustering_plot <- renderPlot({
      hclust_demo <- hclust(dist(matrix(rnorm(100), ncol = 10)))
      plot(hclust_demo, main = "Sample Dendrogram", xlab = "Samples", ylab = "Distance")
    })
    
    output$heatmap_plot <- renderPlot({
      data <- matrix(rnorm(1000), nrow = 100)
      heatmap(data, Colv = NA, Rowv = NA, main = "Methylation Heatmap")
    })
    
    output$global_boxplot <- renderPlot({
      boxplot(list(Sample1 = rnorm(1000, 0.5), Sample2 = rnorm(1000, 0.4), Sample3 = rnorm(1000, 0.6)),
              main = "Mean Methylation by Sample", ylab = "Beta Values")
    })
    
    output$global_density <- renderPlot({
      plot(density(rnorm(1000, 0.5)), main = "Methylation Distribution", xlab = "Beta Values", ylab = "Density")
    })
    
    output$umap_plot <- renderPlot({
      plot(rnorm(50), rnorm(50), main = "UMAP Projection", xlab = "UMAP 1", ylab = "UMAP 2", 
           pch = 16, cex = 2)
    })
    
    # Initialize first view
    show_view("view_beta_matrix", "Beta Matrix")
    update_active_button("nav_beta_matrix")
    
  })
}