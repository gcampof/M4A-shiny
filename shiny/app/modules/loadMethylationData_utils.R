library(limma)

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


# Function to load the raw data and the targets
loadMethylationData <- function(dataDir, arrayType, pattern = "SampleSheet.csv") {
  cat(dataDir)
  
  cat(sprintf("\n Loading %s data ...", arrayType), Sys.time(), "\n")
  dataPath <- file.path(dataDir, arrayType)
  targets <- read.metharray.sheet(dataPath, pattern=pattern)
  rgSet <- read.metharray.exp(dataPath, targets=targets, verbose=TRUE, 
                              force=TRUE, extended=T)
  targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep="_")
  sampleNames(rgSet) <- targets$ID
  index <- seq_along(rgSet@colData@listData)
  rgSet@colData@listData[index] <- lapply(rgSet@colData@listData[index], as.factor)
  cat(sprintf("\n %s data loaded", arrayType), Sys.time(), "\n")
  return(rgSet)
}


reorganizeIdatsIfNeeded <- function(dataDir, all.idats) {
  array_types <- c("450K", "EPIC", "EPIC_V2")
  array_paths <- file.path(dataDir, array_types)
  
  # Create missing directories
  for (i in seq_along(array_paths)) {
    if (!dir.exists(array_paths[i])) {
      dir.create(array_paths[i])
      cat("Created missing folder:", array_paths[i], "\n")
    }
  }
  
  # Categorize IDAT files by size
  idat_sizes <- sapply(all.idats, file.size)
  idat_epic_v2 <- all.idats[idat_sizes > 14000000]
  idat_epic <- all.idats[idat_sizes > 10000000 & idat_sizes <= 14000000]
  idat_450k <- all.idats[idat_sizes <= 10000000]
  
  # Move helper function
  move_files <- function(files, dest_dir) {
    for (file in files) {
      file_name <- basename(file)
      new_path <- file.path(dest_dir, file_name)
      if (!file.exists(new_path)) {
        file.rename(file, new_path)
        cat("Moved", file_name, "to", dest_dir, "\n")
      } else {
        cat("File", file_name, "already exists in", dest_dir, "\n")
      }
    }
  }
  
  # Move files to appropriate folders
  move_files(idat_450k, file.path(dataDir, "450K"))
  move_files(idat_epic, file.path(dataDir, "EPIC"))
  move_files(idat_epic_v2, file.path(dataDir, "EPIC_V2"))
  
  sample_sheet_candidates <- list.files(dataDir, pattern = 'samplesheet\\.csv$', ignore.case = TRUE, full.names = TRUE)
  
  if (length(sample_sheet_candidates) > 0) {
    sample_sheet_path <- sample_sheet_candidates[1]
    cat("\nSampleSheet file found:", basename(sample_sheet_path), "\n")
    samplesheet <- read.csv(sample_sheet_path, stringsAsFactors = FALSE)
    
    # Normalize Array_Type
    if (!"Array_Type" %in% colnames(samplesheet)) {
      cat("ERROR: 'Array_Type' column not found in SampleSheet.csv.\n")
      return()
    }
    
    samplesheet$Array_Type <- toupper(samplesheet$Array_Type)
    samplesheet$Array_Type[samplesheet$Array_Type == "EPIC_V2" | samplesheet$Array_Type == "EPIC_V2"] <- "EPIC_V2"
    
    for (atype in array_types) {
      subset_sheet <- subset(samplesheet, Array_Type == atype)
      out_path <- file.path(dataDir, atype, "SampleSheet.csv")
      if (nrow(subset_sheet) > 0) {
        write.csv(subset_sheet, out_path, row.names = FALSE)
        cat("Wrote", nrow(subset_sheet), "rows to", out_path, "\n")
      } else {
        cat("No rows for", atype, "in SampleSheet\n")
      }
    }
  } else {
    cat("\nWarning: No SampleSheet.csv found in", dataDir, "\n")
  }
  cat("Reorganization complete at", Sys.time(), "\n")
}



normalizeMeth <- function(rgSet) {
  method <- tolower(method)
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




filterDetectionP <- function(rgSet, mSetSq, threshold = detP.threshold) {
  detP <- detectionP(rgSet)
  detP <- detP[match(featureNames(mSetSq), rownames(detP)), ]
  keep <- rowSums(detP < threshold) == ncol(rgSet)
  
  message("Detection p-value filter: ", threshold)
  message("Probes failing: ", sum(!keep), " | Probes kept: ", sum(keep))
  
  mSetSq[keep, ]
}



adjustFFPE <- function(methy, unmethy, tissue_type) {
  unique_types <- unique(tissue_type)
  
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

methylationQC <- function(array_name,rgSet, col_vector, detP.threshold) {
  cat("\nIdentifying failed positions in the methylated and unmethylated channel for the array ", array_name, " at ", Sys.time(), "\n")
  detP <- detectionP(rgSet)
  
  generateBarplot <- function(detP, col_vector, ylim = NULL) {
    barplot(colMeans(detP), col = col_vector[factor(rgSet$Sample_Group)], las = 2,
            cex.names = 0.5, ylab = "Mean detection p-values", names = rgSet$Sample_Name,
            ylim = ylim)
    abline(h = 0.01, col = "red")
    legend("topleft", legend = levels(factor(rgSet$Sample_Group)), fill = col_vector,
           bg = "white", cex = 0.75, ncol = 2)
  }
  
  saveBarplotPDF <- function(filename, detP, col_vector, width, height, ylim = NULL) {
    pdf(filename, width = width, height = height)
    generateBarplot(detP, col_vector, ylim)
    dev.off()
  }
  
  # Generate and save barplots
  print(resultsDir)
  results_file_base <- file.path(paste0(resultsDir,"/1.0-First QC approach for array ", array_name ,".pdf"))
  nb.cols <- length(unique(rgSet$Sample_Group))
  
  saveBarplotPDF(results_file_base, detP, col_vector, width = 18, height = 7)
  saveBarplotPDF(file.path(paste0(resultsDir,"/1.2-BarplotQC_Mean Detection Pvals for array ",array_name,".pdf")), detP, col_vector, width = 44, height = 10)
  
  # Generate QC report
  qcReport(rgSet, sampNames = rgSet$Sample_Name, sampGroups = rgSet$Sample_Group, pdf = file.path(paste0(resultsDir, "/2.0-QC Report",array_name,".pdf")))
  cat("\nQC report created. ", Sys.time(), "\n")
  
  # Check samples with different p-value thresholds
  pval_thresholds <- c(0.1, 0.05, 0.01)
  for (threshold in pval_thresholds) {
    cat(sprintf("\nSamples with a detection p-value (>%.2f): ", threshold), table(colMeans(detP) < threshold), "\n")
    cat("Samples are: ", colnames(detP)[!(colMeans(detP) < threshold)], "\n")
  }
  
  ## ── choose per-sample detection-p threshold ──────────────────────────────
  
  # TODO: ESTO TENGO QUE HACERLO INTERACTIVO SI NO NADA
  cat("\nSelect the sample-level detection-p threshold (proportion of failed probes):\n",
      "  1) 0.10 \n",
      "  2) 0.05 \n",
      "  3) 0.01  \n",
      sep = "")
  ans <- trimws(readline("Enter 1-3 or type the value [0.05]: "))
  
  mapping <- c("1" = 0.10, "2" = 0.05, "3" = 0.01)
  selected_threshold <- if (ans == "") {
    0.05
  } else if (ans %in% names(mapping)) {
    mapping[ans]
  } else if (ans %in% c("0.1", "0.10", "0.05", "0.01")) {
    as.numeric(ans)
  } else {
    cat("Unrecognised input – defaulting to 0.05\n")
    0.05
  }
  cat(sprintf("\nYou have selected a sample detection-p threshold of %.2f\n",
              selected_threshold))
  
  ## ── flag and export per-sample failure rates ─────────────────────────────
  keep <- colMeans(detP) < selected_threshold           # samples to retain
  failed_perc <- colSums(detP > selected_threshold) /
    nrow(detP) * 100                       # % probes failing
  
  cat("\nExporting table with percentage of failed probes per sample …\n")
  outfile <- file.path(
    resultsDir,
    sprintf("1.1-Percentage_of_failed_probes_by_sample_detection_p_%.2f_%s.csv",
            selected_threshold, array_name)
  )
  write.csv(failed_perc,
            file      = outfile,
            quote     = FALSE,
            row.names = TRUE)
  # Exclude poor quality samples
  rgSet <- rgSet[, keep]
  cat(sprintf("\nDimensions of the filtered rgSet: %s\n", dim(rgSet)))
  return(rgSet)
}

plotPostQC <- function(mSet, array_name, resultsDir) {
  pdf(file.path(resultsDir, paste0("3- Post-filtered data QC for array_", array_name, ".pdf")),
      width = 12, height = 7)
  qc <- getQC(mSet)
  plotQC(qc)
  dev.off()
}



processArray <- function(array_name, dataDir=dataDir, filterDir=filterDir, detP.threshold = detP_thr , 
                         normalization.method = norm_method,
                         resultsDir = resultsDir) {
  
  message("Starting array processing: ", array_name)
  message("Detection P-value threshold: ", detP.threshold)
  #message("Normalization method: ", normalization.method)
  
  # Load and QC
  rgSet <- loadMethylationData(dataDir, array_name)

  
  rgSet <- methylationQC(array_name, rgSet, col_vector, detP.threshold)
  
  # Normalization
  mSetSq <- normalizeMeth(rgSet, method =normalization.method)
  
  # Calculate raw & normalized beta values
  betaList <- calculateBeta(rgSet, mSetSq)
  #beta_unorm <- betaList$beta_unorm
  beta_unf <- betaList$beta
  
  # Detection p-value filtering
  # TODO: detp.threshold ES EL THRESHOLD DEL OTRO
  mSetSq <- filterDetectionP(rgSet, mSetSq, threshold = detP.threshold)
  
  # Filter probes
  mSetSqFlt <- filterProbes(mSetSq, codeDir=codeDir)
  
  # Batch correction for FFPE/frozen
  message("Performing FFPE/Frozen correction...")
  beta <- adjustFFPE(
    getMeth(mSetSqFlt),
    getUnmeth(mSetSqFlt),
    mSetSqFlt$Tissue_Type
  )
  
  # Only aggregate to probes for EPIC_V2
  if (array_name == "EPIC_V2") {
    message("Aggregating to probes (EPIC_V2 only)...")
    beta <- aggregate_to_probes(beta)
    #beta_unf <- aggregate_to_probes(beta_unf)
    #beta_unorm <- aggregate_to_probes(beta_unorm)
    mVals_unf <- getM(mSetSq)    
    mVals_unf <- aggregate_to_probes(mVals_unf)
    unfiltered_data <- mSetSq    
    filtered_dat <- mSetSqFlt
  } else {  mVals_unf <- getM(mSetSq)  }
  
  # Final SNP/XY/cross-hyb filtering
  beta <- finalizeBeta(beta)
  
  # QC plot
  plotPostQC(mSetSqFlt, array_name, resultsDir)
  message("Array processing complete: ", array_name)
  saveRDS(beta, file = file.path(resultsDir, paste0("001_beta_", array_name, ".rds")))
  write.csv(beta, file = file.path(resultsDir, paste0("001_beta_", array_name, ".csv")))
  saveRDS(mSetSq, file = file.path(resultsDir, paste0("002_unfilteredData_", array_name, ".rds")))
  saveRDS(mSetSqFlt, file = file.path(resultsDir, paste0("003_filteredData_", array_name, ".rds")))
  
  return(list(
    beta = beta,
    mVals_unf = mVals_unf
    #filtered_data = mSetSqFlt,
    #unfiltered_data = mSetSq,
  ))
}


loadMethylation <- function(dataDir, resultsDir, filterDir, cores, norm_method) {
  
  # ---- basic setup ----
  # Pre-check if legal number of cores
  library(doParallel)
  registerDoParallel(cores)
  
  # # ---- locate IDAT files / array folders ----
  array_types <- c("450K", "EPIC", "EPIC_V2")
  array_dirs  <- setNames(file.path(dataDir, array_types), array_types)
  hasArray    <- sapply(array_dirs, dir.exists)

  available_arrays <- names(hasArray)[hasArray]
  if (length(available_arrays) == 0)
    stop("No valid array folders found after reorganization.")

  # ---- process each array ----
  processed_results <- vector("list", length(available_arrays))
  names(processed_results) <- available_arrays
  
  for (array_type in available_arrays) {
    print('checkpoint during process')
    processed_results[[array_type]] <- processArray(
      array_name           = array_type,
      dataDir              = dataDir,
      filterDir            = filterDir,      # <- now passes through
      detP.threshold       = detP_thr,
      normalization.method = norm_method,
      resultsDir           = resultsDir
    )
  }

  # ---- merge SampleSheets ----
  all_targets <- lapply(available_arrays, function(array_type) {
    dataPath <- file.path(dataDir, array_type)
    targets  <- read.metharray.sheet(dataPath, pattern = "SampleSheet.csv")
    targets$ID <- paste(targets$Sample_Group, targets$Sample_Name, sep = "_")
    targets
  })
  targets <- do.call(rbind, all_targets)
  write.csv(targets, file = file.path(resultsDir, paste0("targets_merged.csv")))
  cat("\nMerged SampleSheet completed at", Sys.time(), "\n")

  # ---- beta-distribution plots ----
  #plot_input <- lapply(processed_results, `[`, c("beta", "beta_unf"))
  #plotBetaValues(plot_input)
  #cat("\nBeta normalization boxplots completed at", Sys.time(), "\n")
  
  print('checkpoint after process')
  
  # ---- merge matrices ----
  merge_matrix <- function(x, y) {
    m <- merge(x, y, by = "row.names")
    rownames(m) <- m[, 1]; m[, 1] <- NULL; m
  }
  
  cat("\nMerging matrices at", Sys.time(), "\n")
  #beta_unf   <- Reduce(merge_matrix, lapply(processed_results, `[[`, "beta_unf"))
  #write.csv(beta_unf, file = file.path(resultsDir, paste0("beta_unfiltered_merged.csv")))
  mVals_unf  <- Reduce(merge_matrix, lapply(processed_results, `[[`, "mVals"))
  #write.csv(mVals_unf, file = file.path(resultsDir, paste0("mVals_unfiltered_merged.csv")))
  beta_merged <- Reduce(merge_matrix,lapply(processed_results, `[[`, "beta"))
  write.csv(beta_merged, file = file.path(resultsDir, paste0("beta_merged.csv")))
  
  mVals <- mVals_unf[rownames(beta_merged), ]
  write.csv(mVals, file = file.path(resultsDir, paste0("mVals_merged.csv")))
  cat("\nFinal merge completed at", Sys.time(), "\n")
  
  return(list(
    beta             = beta_merged,
    #beta_unf         = beta_unf,
    #mVals            = mVals,
    #mVals_unf        = mVals_unf,
    targets          = targets,
    mVals = mVals)
  )
}

merge_matrix <- function(x, y) {
  m <- merge(x, y, by = "row.names")
  rownames(m) <- m[, 1]; m[, 1] <- NULL; m
}



plotBetaValues <- function(beta_list) {
  pdf(paste0(resultsDir,"/2.1-Raw vs. Normalized boxplots.pdf"), width = 17, height = 11) 
  par(mfrow=c(3,2))
  lapply(seq_along(beta_list), function(i) {
    beta_data <- beta_list[[i]]
    boxplot(as.data.frame(beta_data$beta_unorm), main = paste("Original", names(beta_list)[i]))
    boxplot(as.data.frame(beta_data$beta_unf), main = paste("Normalized", names(beta_list)[i]))
  })
  dev.off()
}