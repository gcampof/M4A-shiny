library(shiny)
library(bslib)
library(shinyjs)
source("modules/load_data.R")
source("modules/global_methylation.R")

options(shiny.reactlog = TRUE)
options(shiny.maxRequestSize = 5 * 1024^3)
options(
  ExperimentHub.cache = "/home/rstudio/.cache/R/ExperimentHub",
  ask = FALSE
)

ui <- navbarPage(
  title = "My Pipeline",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  
  header = tagList(shinyjs::useShinyjs()),
  
  tabPanel("Load Data", load_data_ui("load_data")),
  
  tabPanel("Global Methylation", 
           value = "global",
           global_methylation_ui("global_meth")
  ),
  
  tabPanel("Differential Methylation", 
           value = "diffmeth",
           h3("DIFF placeholder")
  )
)

server <- function(input, output, session) {
  
  ready <- reactiveVal(FALSE)
  
  load_data_server("load_data", notify_ready = ready)
  global_methylation_server("global_meth")
  
  # observe({
  #   if (ready()) {
  #     shinyjs::enable(selector = "a[data-value='global']")
  #     shinyjs::enable(selector = "a[data-value='diffmeth']")
  #   } else {
  #     shinyjs::disable(selector = "a[data-value='global']")
  #     shinyjs::disable(selector = "a[data-value='diffmeth']")
  #   }
  # })
}

shinyApp(ui, server)