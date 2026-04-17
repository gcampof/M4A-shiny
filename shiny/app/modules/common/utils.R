# Set up common directories (shared across all analyses)
setup_common_dirs <- function(cfg) {
  common_dirs <- list(
    cache                = here::here(cfg$cache),
    data                 = here::here(cfg$data),
    filter               = here::here(cfg$filter),
    custom_color_palette = here::here(cfg$custom_color_palette),
    pathways             = here::here(cfg$pathways)
  )
  
  for (d in common_dirs) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  
  message("Common directories ready")
  return(common_dirs)
}

# Set up analysis directory (creates new one with auto-generated ID from session)
setup_analysis_dir <- function(common_dirs, cfg, session) {
  analysis_id <- paste0(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    substr(session$token, 1, 8) 
  )
  
  analysis_dir <- file.path(common_dirs$data, paste0("analysis_", analysis_id))
  
  dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
  
  subdirs <- list()
  for (name in names(cfg$subdirs)) {
    sub_path <- file.path(analysis_dir, cfg$subdirs[[name]])
    dir.create(sub_path, showWarnings = FALSE, recursive = TRUE)
    subdirs[[name]] <- sub_path
  }
  
  message("Analysis directory created: ", analysis_dir)
  
  c(common_dirs, list(analysis = analysis_dir, analysis_id = analysis_id), subdirs)
}


cleanup_old_analysis_dirs <- function(base_dir, max_age_hours = 24) {
  if (!dir.exists(base_dir)) return()
  
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  now <- Sys.time()
  
  for (d in dirs) {
    info <- file.info(d)
    if (is.na(info$mtime)) next
    
    age <- difftime(now, info$mtime, units = "hours")
    
    if (age > max_age_hours) {
      unlink(d, recursive = TRUE, force = TRUE)
      message("[CLEANUP] Removed old analysis dir: ", d)
    }
  }
}


# GO: Biological Process gene sets (SYMBOL)
get_go_bp_gene_sets <- function() {
  go_terms <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = keys(org.Hs.eg.db, keytype = "GO"),
    columns = c("GO", "SYMBOL", "ONTOLOGY"),
    keytype = "GO"
  )
  go_bp_terms <- go_terms[go_terms$ONTOLOGY == "BP", ]
  go_bp_gene_sets <- split(go_bp_terms$SYMBOL, go_bp_terms$GO)
  go_bp_gene_sets <- lapply(go_bp_gene_sets, unique)
  go_bp_gene_sets[!sapply(go_bp_gene_sets, is.null)]
}


# Set up annotation cache
setup_cache <- function(DIRS, cfg) {
  cache_dir <- DIRS$cache
  # Built annotation (expensive)
  annot_path <- file.path(cache_dir, "built_annot.rds")
  
  if (!file.exists(annot_path)) {
    message("Building annotation... (this may take a few minutes)")
    built_annot <- methylation_buildannot(cfg$annotation_pkg)
    saveRDS(built_annot, annot_path)
    message("Annotation saved to: ", built_annot)
  } else {
    message("Loading annotation from cache...")
    built_annot <- readRDS(annot_path)
  }
  
  # Raw annotation 
  raw_annot_path <- file.path(cache_dir, "raw_annot.rds")
  
  if (!file.exists(raw_annot_path)) {
    message("Building raw annotation...")
    raw_annot <- as.data.frame(minfi::getAnnotation(cfg$annotation_pkg))
    saveRDS(raw_annot, raw_annot_path)
    message("Raw annotation saved to: ", raw_annot_path)
  } else {
    message("Loading raw annotation from cache...")
    raw_annot <- readRDS(raw_annot_path)
  }
  
  # GO-BP gene sets
  gobp_path <- file.path(cache_dir, "gene_set_list.rds")
  if (!file.exists(gobp_path)) {
    message("Building GO-BP gene sets...")
    gene_set_list <- get_go_bp_gene_sets()
    saveRDS(gene_set_list, gobp_path)
  } else {
    gene_set_list <- readRDS(gobp_path)
  }
  
  pathways_path <- file.path(cache_dir, "pathways.rds")
  kegg_path <- file.path(DIRS$pathways, cfg$gene_set$kegg)
  hallmark_path <- file.path(DIRS$pathways, cfg$gene_set$hallmark)
  
  if (!file.exists(pathways_path)) {
    message("Building pathways...")
    pathways <- list(
      go_bp     = gene_set_list,
      kegg      = fgsea::gmtPathways(kegg_path),
      hallmarks = fgsea::gmtPathways(hallmark_path)
    )
    saveRDS(pathways, pathways_path)
  } else {
    pathways <- readRDS(pathways_path)
  }
  
  message("Setup complete")
  
  list(
    built_annot = built_annot,
    raw_annot   = raw_annot,
    cache_dir   = cache_dir,
    pathways    = pathways
  )
}


# Supported built-in color palettes
get_built_in_color_palettes <- function(){
  return(
    builtin_palettes <- list(
      "Set1 (Brewer)"      = function(n) RColorBrewer::brewer.pal(n, "Set1"),
      "Dark2 (Brewer)"     = function(n) RColorBrewer::brewer.pal(n, "Dark2"),
      "Paired (Brewer)"    = function(n) RColorBrewer::brewer.pal(n, "Paired"),
      "Set2 (Brewer)"      = function(n) RColorBrewer::brewer.pal(n, "Set2"),
      "viridis"            = function(n) viridis::viridis(n),
      "magma"              = function(n) viridis::magma(n),
      "plasma"             = function(n) viridis::plasma(n),
      "cividis"            = function(n) viridis::cividis(n)
    ))
}


# Load custom palettes from directory
load_custom_palettes <- function(dir) {
  txt_files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  
  if (length(txt_files) == 0) return(list())
  
  palettes <- lapply(txt_files, function(f) {
    colors <- readLines(f, warn = FALSE)
    colors <- trimws(colors[nzchar(colors)])
    palette_name <- tools::file_path_sans_ext(basename(f))
    list(name = palette_name, colors = colors)
  })
  
  # Return as named list of functions 
  named <- setNames(
    lapply(palettes, function(p) {
      force(p)
      function(n) p$colors[seq_len(min(n, length(p$colors)))]
    }),
    paste0(sapply(palettes, `[[`, "name"))
  )
  
  named
}


# Prepare all color palettes
prepare_color_palettes <- function(dir) {
  builtin_palettes <- get_built_in_color_palettes()
  custom_palettes <- load_custom_palettes(dir)
  
  all_palettes <- c(builtin_palettes, custom_palettes)
  
  all_palette_choices <- Filter(
    function(x) length(x) > 0,
    list(
      "Custom"   = names(custom_palettes),
      "Built-in" = names(builtin_palettes)
    )
  )
  
  list(
    all_palettes        = all_palettes,
    all_palette_choices = all_palette_choices
  )
}


# Load and validate a new palette
load_new_palette <- function(file_path, palette_name, palette_dir) {
  tryCatch({
    # Read file safely
    if (!file.exists(file_path)) {
      return(list(success = FALSE, message = "File not found"))
    }
    
    # Read lines
    colors <- readLines(file_path, warn = FALSE, encoding = "UTF-8")
    
    # Clean: trim whitespace and remove empty lines
    colors <- trimws(colors)
    colors <- colors[nzchar(colors)]
    
    # Validate: each line must be a valid hex color
    valid_hex <- function(x) {
      grepl("^#([0-9A-Fa-f]{6}|[0-9A-Fa-f]{3})$", x)
    }
    
    invalid_colors <- colors[!valid_hex(colors)]
    
    if (length(invalid_colors) > 0) {
      msg <- paste(
        "Invalid hex colors found:",
        paste(invalid_colors[1:min(3, length(invalid_colors))], collapse = ", ")
      )
      return(list(success = FALSE, message = msg))
    }
    
    if (length(colors) == 0) {
      return(list(success = FALSE, message = "No valid hex colors found in file"))
    }
    
    # Save to palette directory
    dir.create(palette_dir, showWarnings = FALSE, recursive = TRUE)
    output_file <- file.path(palette_dir, paste0(palette_name, ".txt"))
    writeLines(colors, output_file, useBytes = TRUE)
    
    return(list(
      success = TRUE,
      message = paste0("Loaded ", length(colors), " colors from '", palette_name, "'"),
      file = output_file
    ))
    
  }, error = function(e) {
    list(success = FALSE, message = paste("Error:", e$message))
  })
}

start_logging <- function(out_dir) {
  tryCatch({
    log_file_path <- file.path(out_dir, "logs.txt")
    log_conn <- file(log_file_path, open = "wt")
    if (sink.number(type = "output") > 0) sink(type = "output")
    if (sink.number(type = "message") > 0) sink(type = "message")
    sink(log_conn, type = "output", split = TRUE)
    sink(log_conn, type = "message")
    assign(".log_conn", log_conn, envir = .GlobalEnv)
    options(crayon.enabled = FALSE)
    options(cli.num_colors = 1)
    cat("\n=== LOG STARTED", format(Sys.time()), "===\n\n")
  }, error = function(e) {
    warning("[Logging] Failed to start logging: ", e$message)
  })
}


stop_logging <- function() {
  cat("\n=== LOG STOPPED", format(Sys.time()), "===\n")
  
  # Close sinks in reverse order
  if (sink.number(type = "message") > 0) sink(type = "message")
  if (sink.number(type = "output") > 0) sink(type = "output")
  
  # Close and remove the connection
  if (exists(".log_conn")) {
    close(.log_conn)
    rm(.log_conn, envir = .GlobalEnv)
  }
}


load_heavy_components <- function(session, DIRS, cfg, APP_CACHE) {
  showModal(modalDialog(
    title = div(
      class = "text-center",
      icon("microchip", class = "fa-2x mb-2", style = "color: #2c7fb8;"),
      h4("Loading Analysis Components", class = "mt-2")
    ),
    div(
      class = "text-center",
      div(class = "spinner-border text-primary mb-3",
          role = "status",
          style = "width: 3rem; height: 3rem;"),
      p("This may take a moment...", class = "text-muted"),
      p("Loading packages and cache...", class = "small text-muted")
    ),
    footer = NULL,
    easyClose = FALSE,
    size = "s"
  ))
  
  session$onFlushed(function() {
    tryCatch({
      source("modules/common/all_imports.R", local = TRUE)
      cache <- setup_cache(DIRS, cfg)
      APP_CACHE(cache)
      removeModal()
    }, error = function(e) {
      removeModal()
      showModal(modalDialog(
        title = "Error",
        paste("Failed to load components:", e$message),
        footer = modalButton("OK")
      ))
    })
  }, once = TRUE)
  
  invisible(NULL)  # explicit: this function intentionally returns nothing
}