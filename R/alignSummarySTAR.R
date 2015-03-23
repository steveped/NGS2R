#' @title Extracts the alignment summary from STAR log files.
#'
#' @description Extracts the basic summary for fastq files after running fastqc.
#' 
#' @details This will open the log files as output by the alignment tool STAR and load the information into your workspace.
#' The log files in question are the ones traditionally output with the suffix "Log.final.out"
#' 
#' @param logFiles A \code{vector} of log files. Any files which cannot be found will be ignored.
#' If \code{logPath} is not specified, these must contain path information as part of each filename
#' @param logPath The path to the directory where the log files are located. 
#' If no files are specified, this will default to all files which match the pattern "Log.out"
#' @param pattern A text string found in the name of all required logFiles. Not used if the files are specified by name.
#' 
#' @return a list with components \code{$reads, $splicing, $lengths, $indels} and \code{$performance}.
#' Each of thsi will be a \code{data.frame} with samples as rows and the key information as columns.
#' 
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' 
#' @rdname alignSummarySTAR
#' 
#' @export
alignSummarySTAR <- function(logFiles, logPath, pattern="Log.final.out"){
  
  if (missing(logFiles)){
    if(!file.exists(logPath)) stop("\nCannot find the specified path:\n", logPath)
    else logFiles <- list.files(logPath, pattern=pattern, full.names=TRUE) 
    if (length(logFiles) == 0) stop("\nCannot find any files in ", logPath, " which match the pattern ", pattern, "\n")
  }
  else{
    if (!missing(logPath)){
      logFiles <- file.path(logPath, logFiles)
    }
    fc <- sapply(logFiles, file.exists)
    if (length(which(fc)) == 0) stop("\nCannot find any logFiles. Have you specified the path and filenames correctly?\n")
    else logFiles <- logFiles[fc]
  }    
  
  allLogs <- lapply(logFiles, read.delim, header=FALSE, stringsAsFactors=FALSE)
  samps <- gsub("[\\.\\_]$", "", gsub(pattern, "", basename(logFiles)))
  names(allLogs) <- samps
  
  inReads <- sapply(allLogs, function(x){as.integer(x[grep("Number of input reads", x[,1]),2])})
  uniqMap <- sapply(allLogs, function(x){as.integer(x[grep("Uniquely mapped reads number", x[,1]),2])})
  multMap <- sapply(allLogs, function(x){as.integer(x[grep("Number of reads mapped to multiple loci", x[,1]),2])})
  tooManyMap <- sapply(allLogs, function(x){as.integer(x[grep("Number of reads mapped to too many loci", x[,1]),2])})
  notMap <- inReads - uniqMap - multMap - tooManyMap
  reads <- data.frame(inputReads = inReads, uniquelyMapped = uniqMap, multiMapped = multMap,
                      tooManyMapped = tooManyMap, unMapped = notMap)
  
  spl <- t(sapply(allLogs, function(x){as.integer(x[grep("Number of splices", x[,1]),2])}))
  colnames(spl) <- c("Total", "annotatedInSjdb", "GT_AG","GC_AG", "AT_AC", "nonCanonical")

  lengths <- data.frame(averageInput = sapply(allLogs, function(x){x[grep("Average input read length", x[,1]),2]}),
                       averageMapped = sapply(allLogs, function(x){x[grep("Average mapped length", x[,1]),2]}))

  rpb <- sapply(allLogs, function(x){x[grep("rate per base", x[,1]),2]})
  idl <- sapply(allLogs, function(x){as.numeric(x[grep("average length", x[,1]),2])})
  indels <- data.frame(mismatchRate = rpb[1,], 
                       deletionRate = rpb[2,], averageDeletionLength = idl[1,],
                       insertionRate = rpb[3,], averageInsertionLength = idl[2,])
  
  prog <- sapply(allLogs, function(x){x[1:4, 2]})
  progOut <- data.frame(jobStart = as.POSIXlt(prog[1,], format="%b %d %H:%M:%S"),
                        mappingStart = as.POSIXlt(prog[2,], format="%b %d %H:%M:%S"),
                        mappingFinish = as.POSIXlt(prog[3,], format="%b %d %H:%M:%S"),
                        mappingSpeed  = prog[4,])
                        
  
  return(list(reads=reads, splicing=spl, lengths=lengths, indels=indels, performance =progOut))
  
}
#test <- alignSummarySTAR(logPath = getwd())
