
# Build annotation --------------------------------------------------------
#' @name methylation_buildannotation
#' @description Build annotation data.frame to summarise betas from probes 
#' to genes.
#' @param annot Character string with the data base used for annotation. Set
#' to \code{IlluminaHumanMethylationEPICanno.ilm10b4.hg19} by default
#' @value A data.frame with 6 columns: chr, pos, strand, probe, gene, group.
#' 
methylation_buildannot <- function(annot = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"){
  annotation <- as.data.frame(getAnnotation(annot))
  annotation[annotation == ""] <- NA
  annotation[annotation == " "] <- NA
  
  if (annot == "IlluminaHumanMethylationEPICanno.ilm10b2.hg19" | 
      annot == "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"){
    aux <- annotation[, c(1:4, 22:24)]
    aux.long.grp <- strsplit(aux$UCSC_RefGene_Group, ";")
    aux.long.gene <- strsplit(aux$UCSC_RefGene_Name, ";")
    stopifnot(all(sapply(aux.long.gene, length) ==
                    sapply(aux.long.grp, length)))
    aux.len <- sapply(aux.long.grp, length)
    annotation <- data.frame(chr = rep(aux$chr, aux.len),
                             pos = rep(aux$pos, aux.len),
                             strand = rep(aux$strand, aux.len),
                             Name = rep(aux$Name, aux.len),
                             UCSC_Gene = unlist(aux.long.gene),
                             UCSC_Group = unlist(aux.long.grp),
                             stringsAsFactors = FALSE)
    annotation$Group <- gsub("5'UTR|1stExon|TSS200", "Promoter200", annotation$UCSC_Group)
    annotation$Group <- gsub("TSS1500", "Promoter1500", annotation$Group)
    annotation$Group <- gsub("ExonBnd", "Body", annotation$Group)
    annotation <- annotation[!duplicated(annotation), ]
    
    rm(aux, aux.long.gene, aux.long.grp, aux.len)
    annotation <- annotation[, -6]
    annotation <- annotation[! duplicated(annotation), ]
    annotation <- annotation[order(annotation$UCSC_Gene), ]
    return(annotation)
  } else {
    if (annot == "IlluminaHumanMethylation450kanno.ilmn12.hg19"){
      aux <- annotation[, c(1:4, 24:26)]
      aux.long.grp <- strsplit(aux$UCSC_RefGene_Group, ";")
      aux.long.gene <- strsplit(aux$UCSC_RefGene_Name, ";")
      stopifnot(all(sapply(aux.long.gene, length) ==
                      sapply(aux.long.grp, length)))
      aux.len <- sapply(aux.long.grp, length)
      annotation <- data.frame(chr = rep(aux$chr, aux.len),
                               pos = rep(aux$pos, aux.len),
                               strand = rep(aux$strand, aux.len),
                               Name = rep(aux$Name, aux.len),
                               UCSC_Gene = unlist(aux.long.gene),
                               UCSC_Group = unlist(aux.long.grp),
                               stringsAsFactors = FALSE)
      annotation$Group <- gsub("5'UTR|1stExon|TSS200", "Promoter200", annotation$UCSC_Group)
      annotation$Group <- gsub("TSS1500", "Promoter1500", annotation$Group)
      annotation$Group <- gsub("ExonBnd", "Body", annotation$Group)
      annotation <- annotation[! duplicated(annotation), ]
      
      rm(aux, aux.long.gene, aux.long.grp, aux.len)
      annotation <- annotation[, -6]
      annotation <- annotation[! duplicated(annotation), ]
      annotation <- annotation[order(annotation$UCSC_Gene), ]
      return(annotation)
    } else {
      stop("Not implemented yet, only works with IlluminaHumanMethylationEPICanno.ilm10b2/4.hg19 or IlluminaHumanMethylation450kanno.ilmn12.hg19")
    }
  }
}


# Summarise to gene -------------------------------------------------------
#' @name methylation_genemat
#' @description Summarise betas from probes to genes according to position 
#' from TSS.
#' @param beta_mat Beta matrix of probes.
#' @param annot Data frame as returned by function \code{methylation_buildannot}.
#' @param group Grouping factor: TSS200, TSS1500, Body, BodyUTR. TSS200 by 
#' default.
#' @param progressBar Boolean indicating whether or not to print a progress bar.
#' @param rm_mmap Boolean - Whether to remove or not multimapping probes 
#' (i.e. one probe links to more than one group/gene)
#' @value A matrix of beta values per gene and sample.
#' 
methylation_genemat <- function(beta.matrix, annotation, group = "TSS200",
                                progressBar = TRUE, rm_mmap = FALSE){
  group <- match.arg(group, choices = c("TSS200", "TSS1500", "Body", "BodyUTR"))
  
  if (! group %in% c("TSS200", "TSS1500", "Body", "BodyUTR")){
    stop("Group should be one of the stated in the documentation.")
  }
  if (group == "TSS200"){
    groupingfactor <- "Promoter200"
  }
  if (group == "TSS1500"){
    groupingfactor <- c("Promoter200", "Promoter1500")
  }
  if (group == "Body"){
    groupingfactor <- "Body"
  }
  if (group == "BodyUTR"){
    groupingfactor <- c("Body", "3'UTR")
  }
  
  mat <- matrix(0,
                nrow = length(unique(annotation$UCSC_Gene[annotation$Group %in% groupingfactor])),
                ncol = ncol(beta.matrix),
                dimnames = list(c(unique(annotation$UCSC_Gene[annotation$Group %in% groupingfactor])),
                                c(colnames(beta.matrix))))
  if(progressBar){
    pb <- txtProgressBar(min = 0, max = nrow(mat), style = 3)
  }
  if (rm_mmap == TRUE){
    aux <- annotation %>% group_by(Name) %>% 
      summarize(n = length(chr), nd = n_distinct(Group))
    orig_dim <- nrow(aux)
    aux <- aux[aux$nd == 1, ]
    annotation <- annotation[annotation$Name %in% aux$Name, ]
    message("Multimapping probes removed: ", 
            round((100*nrow(aux)/orig_dim), 4),
            "% of probes remaining. ")
    rm(orig_dim)
  } else {
    message("Using all probes")
  }
  for (g in 1:nrow(mat)){
    myProbes <- unique(annotation$Name[which(annotation$UCSC_Gene == rownames(mat)[g] &
                                               annotation$Group %in% groupingfactor)])
    mat[g, ] <- apply(beta.matrix[rownames(beta.matrix) %in% myProbes, , drop = FALSE],
                      2, median, na.rm = TRUE)
    rm(myProbes)
    if (progressBar){
      setTxtProgressBar(pb, g)
    }
  }
  
  mat_clean <- mat[rowSums(is.na(mat)) != (ncol(mat)), ]
  if (progressBar){
    close(pb)
  }
  return(mat_clean)
}

