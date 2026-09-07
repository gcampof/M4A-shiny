# Workers are plain R processes with none of the app loaded, so they set their
# own thread caps and load their own packages. 
m4a_worker_init <- function(app_dir, heavy = TRUE) {
  # The marker lives in globalenv alongside the functions it guards: mirai clears
  # a daemon's globalenv between tasks but R options persist, so an option-based
  # flag would outlive the definitions it vouched for.
  level <- get0(".m4a_worker_ready", envir = globalenv(), ifnotfound = "")
  if (identical(level, "heavy") || (identical(level, "light") && !heavy)) {
    return(invisible(TRUE))
  }

  setwd(app_dir)

  source("modules/common/concurrency.R")
  m4a_apply_thread_caps()

  options(device = function(...) grDevices::pdf(NULL))

  suppressPackageStartupMessages({
    # Attached by app.R in the main process, which a worker never runs.
    library(dplyr)
    library(tibble)
    library(data.table)
    library(matrixStats)

    if (heavy) {
      source("modules/common/all_imports.R")
      library(openxlsx)
      library(ggplot2)
      library(RColorBrewer)
      library(viridis)
      library(colorspace)
    }

    source("modules/common/utils.R")
    source("modules/load_data/load_data_helper.R")

    if (heavy) {
      source("modules/primary_analysis/annotations.R")
      source("modules/primary_analysis/utils.R")
      source("modules/primary_analysis/cnv/cnv_utils.R")
      source("modules/primary_analysis/differential/differential_utils.R")
      source("modules/primary_analysis/heatmap/heatmap_utils.R")
      source("modules/primary_analysis/mds/mds_utils.R")
      source("modules/primary_analysis/pca/pca_utils.R")
      source("modules/primary_analysis/umap/umap_utils.R")
      source("modules/primary_analysis/global_met/global_utils.R")
    }
  })

  assign(".m4a_worker_ready", if (heavy) "heavy" else "light", envir = globalenv())
  gc(full = TRUE)   # loading the stack leaves ~700 MB of reclaimable pages
  invisible(TRUE)
}
