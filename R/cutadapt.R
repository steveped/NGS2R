#' @title A wrapper for the bash shell command cutadapt.
#'
#' @description Takes an individual fastq file then finds it's mate to run as paired end files.
#' 
#' @details This is a specific wrapper which currently gives you no choice about renaming any output files. 
#'  Trimmed Files will be written to the specified directory with identical names to the source files.
#'  As the approach required by \code{cutadapt} for paired ends is to run once on each set of reads, 
#'  this wrapper will write to a folder \code{temp} as a subdirectory of the output directory \code{outDir}.
#'  The temporary files created will be deleted, but the directory will not.
#'  This is to enable easier upscaling to parallel implementations.
#'  
#'  To run in parallel see the examples section
#' 
#'  Note also that the log file contained in the output can be written to the R \code{stdout} in a pretty format using the command \code{cat}.
#'  For example, some important output from the forward set of reads could be printed using \code{cat(log$forward[1:13], sep="\n")}.
#'  
#'  For the cutadapt documentation, please visit \url{http://cutadapt.readthedocs.org/en/latest/}
#'
#' @param inFile the first file in the read pair. Currently needs to contain the pattern as defined by p1 in it's name. The full path is also required
#' @param adapters for \code{cutadapt} this must be a character vector with the adapter sequence(s). 
#' For \code{cutadaptPaired}, this must be a list with two components, forward & reverse, in that order. 
#' Each will be applied to the appropriate read pair.
#' @param outDir the output directory for writing the trimmed fastq file. The names wil be automatically the same as the original file
#' @param type the type of adapter binding. Currently accepts only "-a", "-g" or "-b". Note that the prefix "-" must be included. 
#' For the function \code{cutadaptPaired}, a different type of adapter can be specified for each set of reads.
#' @param args the remaining arguments to the function cutadapt. Defaults to a quality threshold of 20 (-q 20) and a minimum length of 20 (-m 20)
#' @param p1 the pattern specific to the name of the first read inthe pair
#' @param p2 the pattern specific to the name of the second read in the pair
#' @param exec the path to the executable. This may vary from system to system
#' @return If succesfully executed the log file from each pass of cutadapt will be written to the slots \code{forward} and \code{reverse}
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' @seealso \code{\link{cat}}, \code{\link{system}}
#' 
#' 
#' @examples
#' \dontrun{
#' # An example of how to run this in parallel across multiple files
#' adaptSeq <- list(forward = "AGATCGGAAGAGCGGTTC", reverse = "AGATCGGAAGAGCGTCGT")
#' fastqDir <- getwd() # Or wherever the fastq files are
#' outDir <- file.path("trimmed_fastq") # Where you want to write the files
#' R1files <- list.files(fastqDir, pattern="R1") # Get the vector of Read1 files
#'
#' library(snow)
#' cl <- makeSOCKcluster(rep("localhost", 4)) # Make an internal cluster with 4 nodes
#' clusterExport(cl, c("adaptSeq", "outDir")) # export the adapter sequences & the output directory
#' cutadaptLog <- clusterApply(cl, R1files, cutadaptPaired, adapters=adaptSeq, outDir=outDir) # Run in parallel
#' stopCluster(cl)
#' }

#' @rdname cutadapt
#' @export
cutadapt <- function(inFile, adapters, outDir, type="-a", args="-q 20 -m 20", exec="~/.local/bin/cutadapt") {
  
  # Check for the executable
  if(!file.exists(exec)) stop("\nIncorrect location of executable:", exec)
  
  # Check the source file exists
  if (!file.exists(inFile)) stop("File not found: ", inFile)
  
  # Check the structure of the adapters
  if (!is.character(adapters))  stop("adapters must be supplied as a character vector")
  if (length(type) > 1) stop("Type mismatch. You cannot specifiy more than one adapter type.")
  if (!type %in% c("-a", "-g", "-b")) stop("Unknown type. This can only be -a, -g or -b.")
  
  # Create the output directory if required & set the output files
  if (!file.exists(outDir)) dir.create(outDir, recursive=TRUE)
  outFile <- file.path(outDir, basename(inFile))
  
  # Set the adapters to handle multiple sequences
  adapArg <- paste(type, adapters)
  
  # Run the forwards & reverse processes
  cmd <- paste(exec, adapArg, args, "-o", outFile, inFile)
  log <- system(cmd, intern=TRUE)
  
  return(log)
  
}

#' @rdname cutadapt
#' @export
cutadaptPaired <- function(inFile, adapters, outDir, type="-a", args="-q 20 -m 20", p1="R1", p2="R2", exec="~/.local/bin/cutadapt") {
  
  # Check for the executable
  if(!file.exists(exec)) stop("\nIncorrect location of executable:", exec)
  
  # Check the source files exist
  inMate <- gsub(p1, p2, inFile)
  if (!file.exists(inFile)) stop("File not found: ", inFile)
  if (!file.exists(inMate)) stop("File not found: ", inMate)
  
  # Check the structure of the adapters
  if (!is.list(adapters) || length(adapters) < 2) {
    stop("adapters must be supplied as a list with at least two components")
  }
  if(names(adapters)[1]!="forward" || names(adapters)[2]!="reverse") stop("adapters must be supplied as a list with components forward & reverse")
  
  chk <- which(!type %in% c("-a", "-g", "-b"))
  if (length(chk) > 0) stop("Unknown type.", type[chk],"\nThis can only be -a, -g or -b.")
  if (length(type) > 2) stop("Type mismatch. You cannot specify more than one adapter type per direction.")
  type <- rep(type, times=2)
  
  # Create the output directory if required & set the output files
  if (!file.exists(outDir)) dir.create(outDir, recursive=TRUE)
  outFile <- file.path(outDir, basename(inFile))
  outMate <- gsub(p1, p2, outFile)
  
  # Create the temp directory & filenames
  tempDir <- file.path(outDir, "temp")
  if (!file.exists(tempDir)) dir.create(tempDir)
  tempFile <- file.path(tempDir, basename(inFile))
  tempMate <- gsub(p1, p2, tempFile)
  
  # Set the adapters to handle multiple sequences
  typeFwd <- paste(type[1], adapters$forward, collapse=" ")
  typeRev <- paste(type[2], adapters$reverse, collapse=" ")
  
  # Run the forwards & reverse processes
  fwdCmd <- paste(exec, typeFwd, args, "-o", tempFile, "-p", tempMate, inFile, inMate)
  fwd <- system(fwdCmd, intern=TRUE)
  revCmd <- paste(exec, typeRev, args, "-o", outMate, "-p", outFile, tempMate, tempFile)
  rev <- system(revCmd, intern=TRUE)
  
  # Remove the temporary files
  file.remove(c(tempFile, tempMate))
  # Don't remove the temp directory, as this may be required by other threads during 
  # parallel processing
  
  return(list(forward=fwd, reverse=rev))
  
}
