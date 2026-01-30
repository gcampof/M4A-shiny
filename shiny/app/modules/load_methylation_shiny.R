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


# Function to load the raw data and the targets
methylationQC_first_part <- function(input_dir, results_dir, array_name, col_vector) {
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
  
  cat("\nIdentifying failed positions in the methylated and unmethylated channel for the array ", array_name, " at ", Sys.time(), "\n")
  detP <- detectionP(rgSet)
  
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
  return(rgSet)
}


mehtylationQC_second_part <- function(results_dir, rgSet, detP_threshold_input){

  
  # ## ── choose per-sample detection-p threshold ──────────────────────────────
  # 
  # cat("\nSelect the sample-level detection-p threshold (proportion of failed probes):\n",
  #     "  1) 0.10 \n",
  #     "  2) 0.05 \n",
  #     "  3) 0.01  \n",
  #     sep = "")
  # ans <- trimws(readline("Enter 1-3 or type the value [0.05]: "))
  # 
  # mapping <- c("1" = 0.10, "2" = 0.05, "3" = 0.01)
  # selected_threshold <- if (ans == "") {
  #   0.05
  # } else if (ans %in% names(mapping)) {
  #   mapping[ans]
  # } else if (ans %in% c("0.1", "0.10", "0.05", "0.01")) {
  #   as.numeric(ans)
  # } else {
  #   cat("Unrecognised input – defaulting to 0.05\n")
  #   0.05
  # }
  # cat(sprintf("\nYou have selected a sample detection-p threshold of %.2f\n",
  #             selected_threshold))
  
  ## ── flag and export per-sample failure rates ─────────────────────────────
  selected_threshold <- detP_threshold_input
  detP <- detectionP(rgSet)
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


processArray_custom<- function(rgSet, normalization_method, filter_dir){
  # Normalization
  mSetSq <- normalizeMeth(rgSet, method =normalization.method)
  
  # Calculate raw & normalized beta values
  betaList <- calculateBeta(rgSet, mSetSq)
  #beta_unorm <- betaList$beta_unorm
  beta_unf <- betaList$beta
  
  # Detection p-value filtering
  mSetSq <- filterDetectionP(rgSet, mSetSq, threshold = detP.threshold)
  
  # Filter probes
  mSetSqFlt <- filterProbes(mSetSq, codeDir=filter_dir)
  
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
  ))
}