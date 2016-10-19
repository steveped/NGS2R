#' @title Extract Read Information From Cutadapt Log Files
#'
#' @description Takes the output from \code{cutadaptPaired} and extracts key information
#' 
#' @details This will extract the read numbers, the numbers trimmed & the numbers discarded from each file in the paired reads.
#' If the log files are from a set of pairs, it will automatically detect this, and get the information for the entire set of paris
#'
#' @param x The output from \code{cutadaptPaired}
#'
#' @return A \code{data.frame}
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>

#' @rdname cutadaptPairedSummary
#' @export
readCutadaptPairedLogs <- function(x) {
  
  ## Checks need to be added here to make sure it actually is a log file
  
  ## First detect whether is is a single pair or a list of pairs.
  if (length(x)==2 && length(which(c("forward", "reverse") %in% names(x))) == 2) {
    fwdTot <- grep("Processed reads", x$forward, value=TRUE)
    fwdTot <- gsub(".+ ([0-9]+)$", "\\1", fwdTot)
    fwdTrim <- grep("Trimmed reads", x$forward, value=TRUE)
    fwdTrim <- gsub(".+ ([0-9]+) (.+)", "\\1", fwdTrim)
    fwdDisc <- grep("Too short reads", x$forward, value=TRUE)
    fwdDisc <- gsub(".+ ([0-9]+) (.+)", "\\1", fwdDisc)
    revTot <- grep("Processed reads", x$reverse, value=TRUE)
    revTot <- gsub(".+ ([0-9]+)$", "\\1", revTot)
    revTrim <- grep("Trimmed reads", x$reverse, value=TRUE)
    revTrim <- gsub(".+ ([0-9]+) (.+)", "\\1", revTrim)
    revDisc <- grep("Too short reads", x$reverse, value=TRUE)
    revDisc <- gsub(".+ ([0-9]+) (.+)", "\\1", revDisc) 
    return(data.frame(R1Tot = as.integer(fwdTot),
                      R1Trim = as.integer(fwdTrim),
                      R1Discard = as.integer(fwdDisc),
                      R2Tot = as.integer(revTot), 
                      R2Trim = as.integer(revTrim),
                      R2Discard = as.integer(revDisc),
                      Total = as.integer(revTot) - as.integer(revDisc)))
  }
  else { # If it is a list of pairs
    fwdTot <- sapply(x, FUN=function(x){grep("Processed reads", x$forward, value=TRUE)})
    fwdTot <- gsub(".+ ([0-9]+)$", "\\1", fwdTot)
    fwdTrim <- sapply(x, FUN=function(x){grep("Trimmed reads", x$forward, value=TRUE)})
    fwdTrim <- gsub(".+ ([0-9]+) (.+)", "\\1", fwdTrim)
    fwdDisc <- sapply(x, FUN=function(x){grep("Too short reads", x$forward, value=TRUE)})
    fwdDisc <- gsub(".+ ([0-9]+) (.+)", "\\1", fwdDisc)
    revTot <- sapply(x, FUN=function(x){grep("Processed reads", x$reverse, value=TRUE)})
    revTot <- gsub(".+ ([0-9]+)$", "\\1", revTot)
    revTrim <- sapply(x, FUN=function(x){grep("Trimmed reads", x$reverse, value=TRUE)})
    revTrim <- gsub(".+ ([0-9]+) (.+)", "\\1", revTrim)
    revDisc <- sapply(x, FUN=function(x){grep("Too short reads", x$reverse, value=TRUE)})
    revDisc <- gsub(".+ ([0-9]+) (.+)", "\\1", revDisc)
    return(data_frame(Sample = names(fwdTot),
                      R1Tot = as.integer(fwdTot),
                      R1Trim = as.integer(fwdTrim),
                      R1Discard = as.integer(fwdDisc),
                      R2Tot = as.integer(revTot), 
                      R2Trim = as.integer(revTrim),
                      R2Discard = as.integer(revDisc),
                      Total = R2Tot - R2Discard))
  }
  

  
  
}