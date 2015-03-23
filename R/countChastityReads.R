#' @title Find how many reads in a file pass or fail the Illumina Chastity Filter.
#'
#' @description Takes an individual fastq file and counts the reads with the Y/N Chastity flag set
#' 
#' @details Returns a data.frame with four components:
#' 
#' \code{Filename} the name of the fastq file being queried. 
#' 
#' \code{Total} contains the total number of reads in the fastq file.
#' 
#' \code{Pass} contains the number of reads which pass the Illumina filter.
#' 
#' \code{Fail} contains the number of reads which fail the Illumina filter.
#' 
#' @param inFile the name of the file to be queried. Must contain the complete path
#' @return A \code{data.frame} as described above
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' 
#' @rdname countChastityReads

#' @export
countChastityReads <- function(inFile) {
  
  if (length(inFile) > 1) {
    warning("Multiple filenames detected. Only the first file will be queried", immediate.=TRUE)
    inFile <- inFile[1]
  }
  if (!file.exists(inFile)) stop("Could not find the specified file:\n", inFile)
  
  # Detect the file extension
  bs <- basename(inFile)
  ext <- gsub(".+\\.(.+)$", "\\1", bs)
  
  if (ext == "gz") {
    
    cmdTot <- paste("zcat", inFile, "| sed -n 1~4p | wc -l")
    total <- as.integer(system(cmdTot, intern=TRUE))
    cmdPass <- paste("zcat", inFile, "| sed -n 1~4p | egrep -c ' [1-2]:N:[0-9]+:'")
    pass <- as.integer(system(cmdPass, intern=TRUE))
    cmdFail <- paste("zcat", inFile, "| sed -n 1~4p | egrep -c ' [1-2]:Y:[0-9]+:'")
    fail <- as.integer(system(cmdFail, intern=TRUE))
    
  }
  else {
    cmdTot <- paste("sed -n 1~4p", inFile, "| wc -l")
    total <- as.integer(system(cmdTot, intern=TRUE))
    cmdPass <- paste("sed -n 1~4p", inFile, "| egrep -c ' [1-2]:N:[0-9]+:'")
    pass <- as.integer(system(cmdPass, intern=TRUE))
    cmdFail <- paste("sed -n 1~4p", inFile, "| egrep -c ' [1-2]:Y:[0-9]+:'")
    fail <- as.integer(system(cmdFail, intern=TRUE))
  }
  
  if (total != pass + fail) stop("Please check the file manually as the numbers from the default query are incompatible")
  
  out <- data.frame(Filename = bs,
                    Total = total,
                    Pass = pass,
                    Fail = fail)
  return(out)
  
}  
  