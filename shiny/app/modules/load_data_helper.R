# === load_data_helper.R ===
library(minfi)
library(DMRcate)
library(readxl)
library(tools)
library(readr)
library(limma)
library(doParallel)
library(foreach)
library(ggplot2)
library(plotly)
library(DT)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)


ARRAY_SUPPORTED <- list(
  `450K` = c("450k", "450", "hm450", "hm_450", "illumina 450k"),
  EPIC = c("epic","epicv1", "epic_v1", "epic_1", "epic v1", "epicv.1"),
  EPIC_V2 = c("epicv2", "epic_v2", "epic_2", "epic v2", "epicv.2")
)

col_vector<-c(
  "#0d570b", "#54c40a", "#E41A1C", "#d97009", "#52367d", "#874c23", "#bfa21d",
  "#6998b5", "#8c8b8b", "#03465e",  "#e080c3", "#9d83d6", "#14b89c", "#2a4880",
  "#86f793", "#A6761D", "#E31A1C", "#FCCDE5", "#E6AB02", "#F4CAE4", "#FFF2AE",
  "#F781BF", "#CCEBC5", "#8DA0CB", "#E78AC3", "#A6CEE3", "#FFFFCC", "#7570B3",
  "#666666", "#984EA3", "#7FC97F", "#FC8D62", "#CCCCCC", "#CAB2D6", "#F2F2F2",
  "#B2DF8A", "#FDC086", "#FFFF33", "#CCEBC5", "#80B1D3", "#D9D9D9", "#FBB4AE",
  "#4DAF4A", "#66C2A5", "#E6F5C9", "#8DD3C7", "#1F78B4", "#B15928", "#FED9A6",
  "#E7298A", "#D95F02", "#FDDAEC", "#B3B3B3", "#FF7F00", "#FFFFB3", "#A6D854",
  "#33A02C", "#FDB462", "#386CB0", "#BEBADA", "#E5C494", "#B3E2CD", "#FDBF6F",
  "#E5D8BD", "#6A3D9A", "#1B9E77", "#FFED6F", "#B3CDE3", "#A65628", "#FDCDAC",
  "#BEAED4", "#999999", "#F1E2CC"
)



# Creates the new directory and returns the path
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

get_array_types <- function(){
  names(ARRAY_SUPPORTED)
}


get_array_synonims <- function(){
  ARRAY_SUPPORTED
}


is_samplesheet <- function(cols) {
  cols <- tolower(cols)
  any(
    c(
      "array", "slide", "sentrixid", "sentrix_id",
      "sentrixposition", "sentrix_position",
      "sentrixidposition"
    ) %in% cols
  )
}


normalize_array_type <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_character_)
  x <- tolower(trimws(x))
  
  for (array_type in names(ARRAY_SUPPORTED)) {
    # Strict equality against the synonyms
    synonyms <- tolower(trimws(ARRAY_SUPPORTED[[array_type]]))
    if (x %in% synonyms) {
      return(array_type)
    }
  }
  
  NA_character_
}


idat_exists <- function(slide, array, array_type, preprocessing_dir) {
  dir <- file.path(preprocessing_dir, array_type)
  if (!dir.exists(dir)) return(FALSE)
  
  pattern <- paste0("^", slide, "_", array, "_.*\\.idat$")
  
  length(
    list.files(
      dir,
      pattern = pattern,
      ignore.case = TRUE
    )
  ) > 0
}


saveDetectionPBarplotPDF <- function(filename, rgSet, detP, col_vector, width, height, ylim = NULL){
  pdf(filename, width = width, height = height)
  barplot(colMeans(detP), col = col_vector[factor(rgSet$Sample_Group)], las = 2,
          cex.names = 0.5, ylab = "Mean detection p-values", names = rgSet$Sample_Name,
          ylim = ylim)
  abline(h = 0.01, col = "red")
  legend("topleft", legend = levels(factor(rgSet$Sample_Group)), fill = col_vector,
         bg = "white", cex = 0.75, ncol = 2)
  dev.off()
}


parse_idat_files<- function(input_dir, preprocessing_dir) {
  # Find all IDAT files (including subdirectories)
  idat_files <- list.files(input_dir, pattern = "\\.idat$", full.names = TRUE, 
                           recursive = TRUE, ignore.case = TRUE)
  if (length(idat_files) == 0) stop("No IDAT files found in the input directory.")
  
  # Create idats folder
  idats_dir <- create_dir(file.path(preprocessing_dir, "idats"))
  
  # Copy IDATs to idats folder
  for (idat_file in idat_files) {
    file.copy(idat_file, file.path(idats_dir, basename(idat_file)), 
              overwrite = TRUE)
  }
  return(idats_dir)
}


order_idat_per_array <- function(idats_dir, preprocessing_dir){
  if (!dir.exists(idats_dir)) {
    stop("idats_dir does not exist: ", idats_dir)
  }
  
  # Supported arrays
  array_types <- get_array_types()
  array_paths <- file.path(preprocessing_dir, array_types)
  
  # Create missing directories
  for (path in array_paths) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Update idat_files to point to new location
  idat_files <- list.files(idats_dir, pattern = "\\.idat$", full.names = TRUE, 
                           ignore.case = TRUE)
  
  # Size thresholds (bytes)
  SIZE_450K_MAX  <- 10000000
  SIZE_EPIC_MAX  <- 14000000
  
  # Classify by size
  idat_sizes <- sapply(idat_files, file.size)
  idat_450k <- idat_files[idat_sizes <= SIZE_450K_MAX]
  idat_epic <- idat_files[idat_sizes > SIZE_450K_MAX & idat_sizes <= SIZE_EPIC_MAX]
  idat_epic_v2 <- idat_files[idat_sizes > SIZE_EPIC_MAX]
  
  # TODO: This requires solid fallback in case the copy goes wrong needs to be notified
  file.rename(idat_450k, file.path(preprocessing_dir, "450K", basename(idat_450k)))
  file.rename(idat_epic, file.path(preprocessing_dir, "EPIC", basename(idat_epic)))
  file.rename(idat_epic_v2, file.path(preprocessing_dir, "EPIC_V2", basename(idat_epic_v2)))
}


generate_idat_dataframe <- function(preprocessing_dir) {
  array_types <- get_array_types()
  array_paths <- file.path(preprocessing_dir, array_types)
  
  rows <- vector("list", length(array_types))
  
  for (i in seq_along(array_types)) {
    # Find all IDATS in the array folders
    idats <- list.files(
      array_paths[i],
      pattern = "\\.idat$",
      full.names = FALSE,
      ignore.case = TRUE
    )
    
    if (length(idats) == 0) next
    
    prefixes <- unique(
      sub("_(Grn|Red|Gr|Rd)\\.idat$", "", idats, ignore.case = TRUE)
    )
    
    # Extract the slide and the array part
    parts <- strsplit(prefixes, "_")
    
    # Prepare new rows with array information for summary
    rows[[i]] <- data.frame(
      Slide = vapply(parts, `[`, character(1), 1),
      Array = vapply(parts, `[`, character(1), 2),
      Array_Type = array_types[i],
      stringsAsFactors = FALSE
    )
  }
  
  dplyr::bind_rows(rows)
}

# Remove unselected idats to a different folder
separate_unselected_idats <- function(idat_df, selected_indices, preprocessing_dir) {
  
  if (nrow(idat_df) == 0) {
    stop("Empty IDAT dataframe provided.")
  }
  
  if (is.null(selected_indices)) {
    selected_indices <- integer(0)
  }
  
  # Identify unselected samples
  unselected_df <- idat_df[-selected_indices, , drop = FALSE]
  
  if (nrow(unselected_df) == 0) {
    message("All samples selected. Nothing to move.")
    return(invisible(NULL))
  }
  
  # Create not_used folder
  not_used_dir <- create_dir(file.path(preprocessing_dir, "not_used"))
  
  # Loop over unselected samples
  for (i in seq_len(nrow(unselected_df))) {
    
    slide      <- unselected_df$Slide[i]
    array_pos  <- unselected_df$Array[i]
    array_type <- unselected_df$Array_Type[i]
    
    source_dir <- file.path(preprocessing_dir, array_type)
    
    if (!dir.exists(source_dir)) {
      warning("Source directory does not exist: ", source_dir)
      next
    }
    
    # Match both Grn and Red IDATs
    pattern <- paste0("^", slide, "_", array_pos, "_.*\\.idat$")
    
    files_to_move <- list.files(
      source_dir,
      pattern = pattern,
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    if (length(files_to_move) == 0) {
      warning("No IDAT files found for sample: ", slide, "_", array_pos)
      next
    }
    
    ok <- file.rename(
      files_to_move,
      file.path(not_used_dir, basename(files_to_move))
    )
    
    if (!all(ok)) {
      warning(
        "Some IDAT files could not be moved for sample: ",
        slide, "_", array_pos
      )
    }
  }
  
  invisible(TRUE)
}


parse_samplesheets <- function(input_dir, preprocessing_dir) {
  
  # Find candidate files
  ss_files <- list.files(
    input_dir,
    pattern = "\\.(csv|xlsx)$",
    full.names = TRUE,
    ignore.case = TRUE,
    recursive = TRUE
  )
  
  if (length(ss_files) == 0) {
    stop("No CSV/XLSX files found in input directory.")
  }
  
  for (file in ss_files) {
    
    # ---- Read file ----
    df <- tryCatch({
      if (grepl("\\.xlsx$", file, ignore.case = TRUE)) {
        readxl::read_excel(file)
      } else {
        readr::read_csv(file, show_col_types = FALSE)
      }
    }, error = function(e) NULL)
    
    if (is.null(df)) next
    
    cols <- colnames(df)
    
    # ---- Check if SampleSheet? ----
    if (!is_samplesheet(cols)) next
    
    cols_l <- tolower(cols)
    
    # Standardize column access
    slide_col <- cols[cols_l %in% c("slide", "sentrixid", "sentrix_id")][1]
    array_col <- cols[cols_l %in% c("array", "sentrixposition", "sentrix_position")][1]
    array_type_col <- cols[cols_l == "array_type"]
    
    # ---- Step 1: determine array type ----
    array_type <- NA_character_
    
    ## A) From array_type column
    if (length(array_type_col) == 1) {
      for (val in unique(df[[array_type_col]])) {
        array_type <- normalize_array_type(val)
        if (!is.na(array_type)) break
      }
    }
    
    ## B) Probe IDAT folders row-by-row
    if (is.na(array_type)) {
      for (i in seq_len(nrow(df))) {
        slide <- df[[slide_col]][i]
        array <- df[[array_col]][i]
        
        for (atype in get_array_types()) {
          if (idat_exists(slide, array, atype, preprocessing_dir)) {
            array_type <- atype
            break
          }
        }
        if (!is.na(array_type)) break
      }
    }
    
    if (is.na(array_type)) {
      warning("Could not determine array type for SampleSheet: ", basename(file))
      next
    }
    
    # ---- Step 2: filter rows without IDATs ----
    keep <- logical(nrow(df))
    
    for (i in seq_len(nrow(df))) {
      keep[i] <- idat_exists(
        df[[slide_col]][i],
        df[[array_col]][i],
        array_type,
        preprocessing_dir
      )
    }
    
    df_clean <- df[keep, , drop = FALSE]
    
    if (nrow(df_clean) == 0) {
      warning("All rows removed from SampleSheet: ", basename(file))
      next
    }
    
    # ---- Step 3: save SampleSheet ----
    out_dir <- file.path(preprocessing_dir, array_type)
    create_dir(out_dir)
    
    readr::write_csv(
      df_clean,
      file.path(out_dir, "SampleSheet.csv")
    )
  }
  
  invisible(TRUE)
}


saveDetectionPBarplotPDF <- function(filename, rgSet, detP, col_vector, width, height, ylim = NULL){
  pdf(filename, width = width, height = height)
  barplot(colMeans(detP), 
          col = col_vector[factor(rgSet$Sample_Group)], 
          las = 2,
          cex.names = 0.5, 
          ylab = "Mean detection p-values",
          names = rgSet$Sample_Name,
          ylim = ylim)
  abline(h = 0.01, col = "red")
  legend("topleft", legend = levels(factor(rgSet$Sample_Group)), fill = col_vector,
         bg = "white", cex = 0.75, ncol = 2)
  dev.off()
}


generateDetectionPBarplot <- function(array, rgSet, detP, col_vector, interactive , threshold){
  # mean detection p-values for each sample
  mean_detP <- colMeans(detP)
  
  # Create a data frame for plotting
  plot_df <- data.frame(
    Sample = rgSet$Sample_Name,
    MeanDetP = mean_detP,
    Group = factor(rgSet$Sample_Group),
    Color = col_vector[factor(rgSet$Sample_Group)]
  )
  
  # Create the ggplot object
  p <- ggplot(plot_df, aes(x = Sample, y = MeanDetP, fill = Group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = col_vector) +
    ylab("Mean detection p-values")  +
    ggtitle(paste0("Mean detection p-values for ", array))
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    theme(legend.title = element_blank()) 
  
  # Add threshold line if interactive
  if (interactive) {
    p <- p + geom_hline(yintercept = threshold, color = "red", linetype = "dashed")
  }
  
  return(p)
}


load_qc_data_for_arrays <- function(cores, preprocessing_dir, qc_dir){
  array_types <- get_array_types()
  qc_results <- list(
    rgsets = list(),
    detections = list(),
    threshold_selected = list()
  )
  
  for (array in array_types) {
    array_dir <- file.path(preprocessing_dir, array)
    
    # Skip if empty
    if (!dir.exists(array_dir) ||
        length(list.files(array_dir, recursive = TRUE)) == 0) {
      message("Skipping empty or missing array directory: ", array)
      next
    }
    array_qc_dir <- create_dir(file.path(qc_dir, array))
    
    targets <- read.metharray.sheet(array_dir, pattern="SampleSheet.csv")
    rgSet <- read.metharray.exp(base=array_dir, targets=targets, verbose=TRUE, 
                                force=TRUE, extended=TRUE)
    # Prepare metadata
    targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep="_")
    sampleNames(rgSet) <- targets$ID
    index <- seq_along(rgSet@colData@listData)
    rgSet@colData@listData[index] <- lapply(rgSet@colData@listData[index], as.factor)
    
    # Calculate detection p-values
    detP <- detectionP(rgSet)
    
    # Generate and save barplots
    nb.cols <- length(unique(rgSet$Sample_Group))
    p <- generateDetectionPBarplot(array, rgSet, detP, col_vector, FALSE , 0)
    barplot_title <- file.path(array_qc_dir, paste0("1.0-BarplotQC_Mean_Detection_Pvals_array_", array, ".pdf"))
    ggsave(barplot_title, p, width = 18, height = 7)
    
    # Generate QC report
    qcReport(rgSet, sampNames = rgSet$Sample_Name, sampGroups = rgSet$Sample_Group,
             pdf = file.path(paste0(array_qc_dir, "/2.0-QC Report",array,".pdf")))
    cat("\nQC report created. ", Sys.time(), "\n")
    
    # Store in results
    qc_results$rgsets[[array]] <- rgSet
    qc_results$detections[[array]] <- detP
    
    # Check samples with different p-value thresholds
    pval_thresholds <- c(0.1, 0.05, 0.01)
    for (threshold in pval_thresholds) {
      cat(sprintf("\nSamples with a detection p-value (>%.2f): ", threshold), table(colMeans(detP) < threshold), "\n")
      cat("Samples are: ", colnames(detP)[!(colMeans(detP) < threshold)], "\n")
    }
  }
  return(qc_results)
}


normalizeMeth <- function(rgSet, norm_method) {
  method <- tolower(norm_method)
  message("Normalizing using method: ", method)
  
  switch(method,
         ssnoob    = preprocessNoob(rgSet),
         raw       = preprocessRaw(rgSet),
         illumina  = preprocessIllumina(rgSet),
         quantile  = preprocessQuantile(rgSet),
         funnorm   = preprocessFunnorm(rgSet),
         stop("Unknown normalization method: ", method,
              "\nValid options are: ssnoob, raw, illumina, quantile, funnorm"))
}


calculateBeta <- function(rgSet, mSetSq) {
  methy_unorm <- getMeth(preprocessRaw(rgSet))
  unmethy_unorm <- getUnmeth(preprocessRaw(rgSet))
  beta_unorm <- methy_unorm / (methy_unorm + unmethy_unorm + 100)
  
  methy <- getMeth(mSetSq)
  unmethy <- getUnmeth(mSetSq)
  beta <- methy / (methy + unmethy + 100)
  
  list(
    beta = beta,
    beta_unorm = beta_unorm,
    methy = methy,
    unmethy = unmethy
  )
}

filterDetectionP <- function(rgSet, mSetSq, threshold) {
  detP <- detectionP(rgSet)
  detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
  keep <- rowSums(detP < threshold) == ncol(rgSet)
  
  message("Detection p-value filter: ", threshold)
  message("Probes failing: ", sum(!keep), " | Probes kept: ", sum(keep))
  
  mSetSq[keep, ]
}


filterProbes <- function(mSetSq, filter_dir) {
  readFilt <- function(fname) read.table(file.path(filter_dir, fname), header = FALSE)[, 1]
  
  amb <- readFilt("amb_3965probes.vh20151030.txt")
  epic <- readFilt("epicV1B2_32260probes.vh20160325.txt")
  snp <- readFilt("snp_7998probes.vh20151030.txt")
  xy <- readFilt("xy_11551probes.vh20151030.txt")
  
  rs <- grep("rs", rownames(mSetSq), value = TRUE)
  ch <- grep("ch", rownames(mSetSq), value = TRUE)
  
  remove <- unique(c(amb, epic, snp, xy, rs, ch))
  keep <- !rownames(mSetSq) %in% remove
  
  message("Filtering probes: removed ", sum(!keep), ", kept ", sum(keep))
  mSetSq[keep, ]
}


adjustFFPE <- function(methy, unmethy, tissue_type) {
  unique_types <- unique(tissue_type)
  
  if (is.null(tissue_type)) {
    message("No FFPE/Frozen adjustment applied (no tissue type column provided in Sample Sheet")
    return(methy / (methy + unmethy + 100))
  }
  
  if (length(unique_types) == 1) {
    message("No FFPE/Frozen adjustment needed (only one tissue type: ", unique_types, ")")
    return(methy / (methy + unmethy + 100))
  }
  
  message("Adjusting by FFPE/Frozen tissue type...")
  batch <- ifelse(tissue_type == "FFPE", 2, 1)
  methy.ba <- 2^removeBatchEffect(log2(methy + 1), batch)
  unmethy.ba <- 2^removeBatchEffect(log2(unmethy + 1), batch)
  
  methy.ba / (methy.ba + unmethy.ba + 100)
}


finalizeBeta <- function(beta) {
  beta <- aggregate_to_probes(beta)
  rmSNPandCH(beta, dist = 2, mafcut = 0.05, rmXY = TRUE, rmcrosshyb = TRUE)
}


plotPostQC <- function(mSet, array_name, resultsDir) {
  pdf(file.path(resultsDir, paste0("3- Post-filtered data QC for array_", array_name, ".pdf")),
      width = 12, height = 7)
  qc <- getQC(mSet)
  plotQC(qc)
  dev.off()
}

# do this for each array_type
generate_beta_matrix <- function(cores, array, rgSet, detP, norm_method, threshold, 
                                 filter_dir, beta_dir) {
  message("Processing array: ", array)
  array_beta_dir <- create_dir(file.path(beta_dir, array))
  
  ## --- 0. Exclude poor quality samples ---
  keep <- colMeans(detP) < threshold
  rgSet <- rgSet[, keep]
  
  ## ---- 1. Normalization ----
  mSetSq <- normalizeMeth(rgSet, norm_method)
  
  ## ---- 2. Raw & normalized beta/M values ----
  beta_list <- calculateBeta(rgSet, mSetSq)
  beta_unf <- beta_list$beta
  
  ## ---- 3. Detection p-value filtering ----
  mSetSq <- filterDetectionP(rgSet, mSetSq, threshold)
  
  ## ---- 4. Probe filtering ----
  mSetSq_flt <- filterProbes(mSetSq, filter_dir)
  
  ## ---- 5. FFPE / Frozen adjustment ----
  beta <- adjustFFPE(getMeth(mSetSq_flt), getUnmeth(mSetSq_flt), mSetSq_flt$Tissue_Type)
  
  ## ---- 6. Array-specific handling ----
  if (array == "EPIC_V2") {
    message("EPIC_V2 detected → aggregating to probes")
    beta <- aggregate_to_probes(beta)
    mVals_unf <- getM(mSetSq)    
    mVals_unf <- aggregate_to_probes(mVals_unf)
    unfiltered_data <- mSetSq    
    filtered_dat <- mSetSq_flt
  } else {
    mVals_unf <- getM(mSetSq)
  }
  
  ## ---- 7. Final SNP / XY / cross-hyb filtering ----
  beta <- finalizeBeta(beta)
  
  ## ---- 8. QC plots ----
  plotPostQC(mSetSq_flt, array, array_beta_dir)
  
  ## ---- 9. Save outputs ----
  saveRDS(beta, file.path(array_beta_dir, paste0("001_beta_", array, ".rds")))
  write.csv(beta, file.path(array_beta_dir, paste0("001_beta_", array, ".csv")))
  saveRDS(mSetSq,      file.path(array_beta_dir, paste0("002_unfilteredData_", array, ".rds")))
  saveRDS(mSetSq_flt,  file.path(array_beta_dir, paste0("003_filteredData_", array, ".rds")))
  
  message("Finished array: ", array)
  
  return(list(
    beta = beta,
    mVals_unf = mVals_unf
  ))
}


merge_samplesheets <- function(arrays, preprocessing_dir, beta_merge_dir){
  all_targets <- lapply(arrays, function(array) {
    array_path <- file.path(preprocessing_dir, array)
    targets <- read.metharray.sheet(array_path, pattern = "SampleSheet.csv")
    targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep = "_")
    targets
  })
  
  # Get all unique columns
  all_cols <- unique(unlist(lapply(all_targets, colnames)))
  
  # Add missing columns as NA to each data frame
  all_targets <- lapply(all_targets, function(df) {
    for (col in all_cols) {
      if (!col %in% colnames(df)) {
        df[[col]] <- NA
      }
    }
    df[, all_cols]
  })
  
  targets <- do.call(rbind, all_targets)
  write.csv(targets, file = file.path(beta_merge_dir, paste0("targets_merged.csv")))
  cat("\nMerged SampleSheet completed at", Sys.time(), "\n")
}

merge_matrix <- function(x, y) {
  m <- merge(x, y, by = "row.names")
  rownames(m) <- m[, 1]; m[, 1] <- NULL; m
}

merge_beta_matrix <- function(processed_results, beta_merge_dir){
  print(processed_results)
  beta_merged <- Reduce(merge_matrix,lapply(processed_results, `[[`, "beta"))
  write.csv(beta_merged, file = file.path(beta_merge_dir, paste0("beta_merged.csv")))
  
  mVals_unf  <- Reduce(merge_matrix, lapply(processed_results, `[[`, "mVals_unf"))
  mVals <- mVals_unf[rownames(beta_merged), ]
  write.csv(mVals, file = file.path(beta_merge_dir, paste0("mVals_merged.csv")))
  cat("\nMerged BetaMatrix completed at", Sys.time(), "\n")
  
}

extract_beta_and_targets <- function(input_dir, beta_dir){
  
  # Get all CSV files recursively from input directory
  all_files <- list.files(path = input_dir, pattern = "\\.csv$", 
                          recursive = TRUE, full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("No CSV files found in the input folder")
  }
  
  beta_matrix_file <- NULL
  targets_file <- NULL
  
  # Check each CSV file to validate if it's a beta matrix or targets file
  for (file_path in all_files) {
    tryCatch({
      df <- read.csv(file_path, nrows = 5)  # Read first few rows for validation
      
      # Check if it's a beta matrix file
      # Typically has numeric columns (except possibly first column with probe names)
      if (is.null(beta_matrix_file)) {
        # Check if it looks like a beta matrix:
        # - Multiple columns (at least 2)
        # - Mostly numeric data
        # - Often has row names/probe IDs
        if (ncol(df) >= 2) {
          numeric_cols <- sapply(df, is.numeric)
          numeric_proportion <- sum(numeric_cols) / ncol(df)
          
          # If at least 80% of columns are numeric, likely a beta matrix
          if (numeric_proportion >= 0.8) {
            beta_matrix_file <- file_path
          }
        }
      }
      
      # Check if it's a targets file
      # Typically has metadata columns like Sample_Name, Basename, Sentrix_ID, etc.
      if (is.null(targets_file)) {
        col_names_lower <- tolower(colnames(df))
        
        # Look for typical targets file columns
        targets_markers <- c("sample_name", "basename", "sentrix", "chip", 
                             "position", "array", "slide")
        
        has_targets_markers <- any(targets_markers %in% col_names_lower)
        
        if (has_targets_markers) {
          targets_file <- file_path
        }
      }
      
    }, error = function(e) {
      # Skip files that can't be read
      NULL
    })
  }
  
  # Validate that both files were found
  if (is.null(beta_matrix_file)) {
    stop("Beta matrix file not found in input directory: ", input_dir)
  }
  
  if (is.null(targets_file)) {
    stop("Targets file not found in input directory: ", input_dir)
  }
  
  # Copy files to output directory
  file.copy(from = beta_matrix_file, 
            to = file.path(beta_dir, "beta.csv"), overwrite = TRUE)
  
  file.copy(from = targets_file, 
            to = file.path(beta_dir, "targets.csv"), overwrite = TRUE)
  
  message("Successfully copied files:")
  message("  Beta matrix: ", basename(beta_matrix_file))
  message("  Targets file: ", basename(targets_file))
  
  invisible(list(
    beta_path = file.path(beta_dir, "beta.csv"),
    targets_path = file.path(beta_dir, "targets.csv")
  ))
}

#TODO: Preguntar a Carla
generate_beta_qc <- function(beta_path, targets_path, beta_dir){
  # Generate QC plots
  beta <- read.csv(beta_path)
  targets <- read.csv(targets_path)
  out_file <- file.path(beta_dir, paste0('M4A_QC_', format(Sys.Date(), '%Y%m%d'), '.pdf'))
  rgset <- NULL
  rgset <- minfi::read.metharray.exp(base = dirname(beta_path))
  
  # TODO: Esto no soy capaz de hacerlo funcionar
  # methylation_QCplots(rgset, file = out_file, sampGroups = targets$Sample_Group, palette = col_vector)
  message("Succesfully generated QC plots for beta matrix")
}