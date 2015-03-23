#' @title Remove reads from a fastq file that fail the Illumina Chastity Filter
#'
#' @description A wrapper to the tool \code{fastq_illumina_filter} which must be installed on the local machine.
#' 
#' @details Returns the log file from the process
#' 
#' @param inFile the name of the file to be processed. Must contain the complete path
#' @param outDir the name of the directory to write the files to.
#' @param args any additional arguments to specify
#' @param exec the path to the executable
#' 
#' @return The text from the logfile, from the command line tool \code{fastq_illumina_filter}
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' 
#' @rdname fastqIlluminaFilter

#' @export
fastqIlluminaFilter <- function(inFile, outDir, args="-N -v", exec="/usr/local/bin/fastq_illumina_filter") {
  
  # Check for the executable
  if(!file.exists(exec)) stop("\nIncorrect location of executable:", exec)
  
  if (length(inFile) > 1) {
    warning("Multiple filenames detected. Only the first file will be queried", immediate.=TRUE)
    inFile <- inFile[1]
  }
  if (!file.exists(inFile)) stop("Could not find the specified file:\n", inFile)
  
  # Detect the file extension
  bs <- basename(inFile)
  ext <- gsub(".+\\.(.+)$", "\\1", bs)
  
  # Check that the output directory exists
  if (!file.exists(outDir)) dir.create(outDir, recursive=TRUE)
  outFile <- file.path(outDir, bs)
  
  if (ext != "gz") cmd <- paste(exec, args, "-o", outFile, inFile)
  else  cmd <- paste("zcat", inFile, "|", exec, args, "| gzip >", outFile)
  
  cat("Filtering reads from", inFile, "\n")
  log <- system(cmd, intern=TRUE)

  return(log)
  
}