#' @title Extract the PASS/FAIL information from a batch of FASTQC summaries
#' 
#' @description Look through a batch of .zip files produced by FASTQC to obtain summary information across a range of files
#' 
#' @details This will look in a number of .zip files, as produced by FASTQC, 
#' and will extract the PASS/FAIL/WARN information across the key report sections.
#' 
#' \code{fqcSummary} will return the PASS/FAIL/WARN information for any number of FASTQC reports.
#' 
#' \code{fqcSummaryPlot} will plot the output from the function \code{fqcSummary} using default parameters.
#' 
#' \code{extractFqcSummary} is the low-level function for a single file and is primarily to be
#' called internally by the other two functions.
#' 
#' @param fqName the filename to extract the totals for. 
#' This should be the name of a single fastq file. 
#' @param fqNames the filenames to extract the totals for.
#' @param qcDir the directory to look in for the FASTQC reports
#' 
#' @return the PASS/FAIL summary information from the FASTQC report(s) either as a 
#' data.frame (\code{fqcSummary}) or as a plot (\code{fqcSummaryPlot})
#' 
#' @import ggplot2
#' @import dplyr
#' @import reshape2
#' 
#' @seealso fastqc
#' 
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @rdname fqcSummary
#' 
extractFqcSummary <- function(fqName, qcDir){
  
  fullName <- list.files(qcDir, fqName, full.names=TRUE)
  zipFile <- grep(".zip$", fullName, value=TRUE)
  if (length(zipFile)>1) stop("The given fqName matches more than one file: ", fqName)
  
  intDir <- gsub(".zip", "", basename(zipFile))
  datFile <- file.path(intDir, "summary.txt")
  data <- readLines(uz <- unz(zipFile, datFile), 12L)
  close(uz)
  
  out.vec <- unlist(strsplit(data, split="\t"))
  out <- out.vec[seq(1, length(out.vec), by=3)]
  names(out) <- out.vec[seq(2, length(out.vec), by=3)]
  
  return(out)
  
}

#' @rdname fqcSummary
#' @export
fqcSummary <- function(fqNames, qcDir){
  
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
  
  counts <-  data.frame(lapply(fqNames, FUN=extractFqcSummary, qcDir=qcDir))
  colnames(counts) <- gsub("_fastqc", "", fqNames)
  return(counts)
  
}

#' @rdname fqcSummary
#' @export
fqcSummaryPlot <- function(fqNames, qcDir){
  
  counts <- fqcSummary(fqNames, qcDir)
  cats <- rownames(counts)
  counts <- mutate(counts, category=factor(cats, levels=cats[seq(length(cats), 1, by=-1)]))
  n <- ncol(counts)
  
  ggCounts <- melt(counts, id.vars="category", measure.vars=colnames(counts)[-n])
  ang <- ifelse(n>3, 90, 0)
  ttl <- paste("FASTQC Summary\n", qcDir, sep="")
  ggplot(ggCounts, aes(x=variable, y=category, fill=value)) +
    geom_tile(colour="black") +
    scale_fill_manual(values=c(FAIL="red",PASS="green",WARN="yellow")) +
    theme(axis.text.x = element_text(angle=ang)) +
    labs(x="Source File", y="QC Category", title=ttl) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0))
  
}