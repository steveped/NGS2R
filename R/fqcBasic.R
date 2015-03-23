#' @title Extracts the basic summary from fastqc reports.
#'
#' @description Extracts the basic summary for fastq files after running fastqc.
#' 
#' @details Always use the function \code{fqcBasic} as this will extract the basic summary information for one or many 
#'  fastqc reports. The function \code{extractFqcBasic} is the low-level function called by the function \code{fqcBasic} 
#'  and is designed to extract information from only a single \code{FASTQC} report.
#' 
#' @param fqName the filename to extract the totals for. 
#' This should be the name of a single fastq file. 
#' @param fqNames a vector of filenames to extract the totals for.
#' @param qcDir the directory to look for the FASTQC reports
#' @return a \code{data.frame} with the basic summary information as contained in FASTQC reports
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @import dplyr

#' @rdname fqcBasic
#' @export
fqcBasic <- function(fqNames, qcDir){
  
  if (missing(qcDir)) stop("The directory must be specified!")
  if (!file.exists(qcDir)) stop("Cannot find the requested directory: ", qcDir)
  if(missing(fqNames)) {
    fqNames <- unique(gsub(".(zip|html)$", "", list.files(qcDir)))
  }
  else {
    allFiles <- list.files(qcDir, pattern=".zip")
    allMatches <- sapply(fqNames, FUN=grep, x=allFiles, value=TRUE)
    if (length(unlist(allMatches)) < length(fqNames)) stop("Some requested files couldn't be found")
    if (length(unlist(allMatches)) > length(fqNames)) stop ("Some requested names matched to multiple files")
  }

  summary <-  lapply(fqNames, FUN=extractFqcBasic, qcDir=qcDir)
  summary <- suppressWarnings(bind_rows(summary))
  return(summary)
  
}

#' @rdname fqcBasic
#' @export
extractFqcBasic <- function(fqName, qcDir){
  
  fullName <- list.files(qcDir, fqName, full.names=TRUE)
  zipFile <- grep(".zip$", fullName, value=TRUE)
  if (length(zipFile)>1) stop("The given fqName matches more than one file: ", fqName)
  
  intDir <- gsub(".zip", "", basename(zipFile))
  datFile <- file.path(intDir, "fastqc_data.txt")
  data <- readLines(uz <- unz(zipFile, datFile), 10L)
  close(uz)
  data <- unlist(strsplit(data[-(1:3)], split="\t"))
  out <- data[seq(2, length(data), by=2)]
  i <- suppressWarnings(sapply(out, as.integer))
  names(out) <- data[seq(1, length(data), by=2)]
  out <- data.frame(t(out))
  out[,!is.na(i)] <- i[!is.na(i)]
  return(out)
  
}