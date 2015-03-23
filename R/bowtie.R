#' @title A wrapper for the bash shell tool bowtie.
#'
#' @description Builds a bowtie index, or aligns to the genome using bowtie.
#' 
#' @details \code{bowtieBuild} will build the index files required for aligning using bowtie, or tophat.
#' The list returned by this function will contain the \code{stdout} information in the component \code{$log}.
#' The \code{ebwt} prefix required by the subsequent functions will be in the component \code{$ebwt}.
#' The remaining list components contain the \code{system} command as executed, the version information and the date the index was built.
#' 
#' Be aware that this can be a lengthy process depending on the genome size. 
#' Please be patient and check the output folder to see if files are still being written to.
#' 
#' A wrapper for bowtie is still under construction and hasn't been written yet
#' 
#' @param genome The reference genome as a single fasta file.
#' @param ebwt The prefix to give the index files. Will default to the basename of the reference genome.
#' @param genDir The directory containing the reference genome. Only required if not specified as part of the genome filename.
#' @param inDir The directory to write the index files to. Will default to the same directory as the reference genome
#' @param args Any additional arguments/options to pass to the commad
#' @param exec the path to the executable. This may vary from system to system
#' @return a \code{list} with components \code{$log}, \code{$ebwt}, \code{$version}, \code{$command} and \code{$time}
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' 
#' @rdname Bowtie
#' @export
buildBowtie <- function(genome, ebwt, genDir, indDir, args=c(), exec="/usr/bin/bowtie-build") {
  
  # Check for the executable
  if(!file.exists(exec)) stop("\nIncorrect location of executable:", exec)
  
  # Check for the reference genome
  if (missing(genDir)) genDir <- dirname(genome)
  genFile <- file.path(genDir, basename(genome))
  if (!file.exists(genFile)) stop("\nCannot locate genome:", genFile)
  if (gsub(".+\\.(fa)$", "\\1", genFile) != "fa") stop("\nGenome must be in fasta format with suffix fa")
  cat("Found reference genome:", genFile, "\n")
  
  # Set the prefix for the ebwt index files
  if (missing(ebwt)) ebwt <- gsub("(.+)(\\..+$)", "\\1", basename(genome))
  cat("Bowtie index will be created using the prefix", ebwt, "\n")
  
  # Set the write location for the index files
  if (missing(indDir)) indDir <- genDir
  if (!file.exists(indDir)) dir.create(indDir, recursive=TRUE)
  indPre <- file.path(indDir, ebwt)
  cat("Index files will be written to the directory", indDir, "\n")

  # Check the additional arguments
  if (!is.null(args)){
    if (length(grep(" ", args)) == 0) {
      argVec <- args
    }
    else {
      argVec <- unlist(strsplit(args, " "))
    }
    argFail <- which(!argVec %in% paste("-", 
                                        c("f", "c", "C", "-color", "a", "-noauto", "p", "-packed",
                                          "B", "-ndoc", "r", "-noref", "3", "-justref", "-ntoa", "q"),
                                        sep=""))
    if (length(argFail) >= 1) {
      stop("You have specified options which are not yet implemented in this wrapper.\n",
           "Please refer to the help page accessed by bowtie-build -h, and run the command from the terminal")                                
    }
    else{
      cat("Specified arguments match those implemented in this wrapper.\n")
    }
  }
  
  st <- Sys.time()
  cat("Commencing index build at", date(), "\n")
  
  vers <- system(paste(exec, "--version"), intern=TRUE)
  cat("Building index using:", vers[1], "\n")
  
  cmd <- paste(exec, args, genFile, indPre)
  log <- system(cmd, intern=TRUE)
  
  ft <- Sys.time()
  return(list(log = log, ebwt = indPre, version = vers, command = cmd, time = data.frame(start=st, finish=ft)))
  
}

# bowtie <- function(ebwt, exec="/usr/bin/bowtie-build") {
#   
# }



