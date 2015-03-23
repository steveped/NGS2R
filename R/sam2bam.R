#' @title Convert SAM files to BAM files
#'
#' @description Converts SAM files to BAM files using \code{samtools view}
#' 
#' @details This is a wrapper to the function \code{samtools view} set to the conversion of sam files to bam files.
#' 
#' Files will be converted and the names of the newly created bam files will be supplied in the output in the component \code{$outFiles}.
#' 
#' This function is designed for simple operations on files existing on a local HDD.
#' For more complex analyses or processes, please see the package \code{Rsamtools} on 
#' Bioconductor (\url{http://bioconductor.org/packages/release/bioc/html/Rsamtools.html}).
#' 
#' @param inFiles A character vector with the full path of the SAM file(s)
#' @param outFiles The names to give the output files. 
#' If not specified, this will default to the same names as the input, but with the suffix .bam
#' @param cl A cluster as defined by the package \code{snow}. If specified, the conversion will be run in parallel
#' @param exec The path to the samtools executable. Defaults to "/usr/bin/samtools". 
#' This may vary depending on where you have installed \code{samtools}
#' 
#' @return a list with the components \code{$inFiles} & \code{$outFiles}
#' 
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' 
#' @rdname sam2bam
#' 
#' @import parallel
#' 
#' @export
sam2bam <- function(inFiles, outFiles, cl=NULL, sort=TRUE, sortArgs=c(), sortPre = "Sorted", rmUnsort=TRUE, index=TRUE, exec= "/usr/bin/samtools"){
  
  fEx <- file.exists(inFiles)
  if (length(which(fEx)) == 0) stop("\nNo files were specified that exist")
  if (length(which(!fEx)) > 0 ) {
    warning("\nCould not find the following file(s), and they will be ignored:\n", inFiles[!fEx] )
    inFiles <- inFiles[fEx]
  }
  
  areSam <- grep("[Ss][Aa][Mm]$", inFiles)
  if (length(areSam) < length(inFiles)) {
    warning("\nSome specified files appear to not be .sam files. These will be ignored\n")
    inFiles <- inFiles[areSam]
  }
  
  if (missing(outFiles)){
    outFiles <- gsub("[Ss][Aa][Mm]$", "bam", inFiles)
  }
  else {
    if (length(outFiles) != length(inFiles)) stop("The supplied outFiles vector does not match the list of input files.\n")
    notBam <- grep("[Bb][Aa][Mm]$", outFiles, invert=TRUE)
    outFiles[notBam] <- paste(outFiles[notBam], ".bam", sep="")
  }
  
  nF <- length(outFiles)
  cmd <- paste(exec, "view -bS -o", outFiles, inFiles)

  if (is.null(cl)) log <- sapply(cmd, FUN=system, intern=TRUE, USE.NAMES=FALSE)
  else log <- parSapply(cl, cmd, system, intern=TRUE, USE.NAMES=FALSE)
  fails <- which(unlist(log)==1)
  if (length(fails)>0) warning("\nThe conversion has failed on:\n", inFiles[fails])
  
  if (sort){ # If sorted files are required
    sortFiles <- gsub("bam", sortPre, outFiles)
    cmd <- paste(exec, "sort", sortArgs, outFiles, sortFiles)
    if (is.null(cl)) log <- sapply(cmd, FUN=system, intern=TRUE, USE.NAMES=FALSE)
    else log <- parSapply(cl, cmd, system, intern=TRUE, USE.NAMES=FALSE)
    fails <- which(unlist(log)==1)
    if (length(fails)>0) warning("\nThe sorting has failed on:\n", outFiles[fails])
    if (rmUnsort) {
      file.remove(outFiles)
      outFiles <- c()
    }
  }
  
  if (index){
    cmd <- paste(exec, "index", paste(sortFiles, "bam", sep="."))
    if (is.null(cl)) log <- sapply(cmd, FUN=system, intern=TRUE, USE.NAMES=FALSE)
    else log <- parSapply(cl, cmd, system, intern=TRUE, USE.NAMES=FALSE)
    fails <- which(unlist(log)==1)
    if (length(fails)>0) warning("\nThe indexing has failed on:\n", sortFiles[fails])
    indexFiles <- paste(sortFiles, "bam", "bai", sep=".")
    fEx <- file.exists(indexFiles)
    if (length(which(!fEx)) > 0) warning("\nError - The Index file cannot be found for:\n", indexFiles[!fEx])
    indexFiles <- indexFiles[fEx]
  }
  
  return(list(inFiles = inFiles, unsorted = outFiles, sorted=paste(sortFiles, "bam", sep="."), indexFiles=indexFiles))
  
}
# cl <- makeCluster(8, type = "SOCK")
# test <- list.files(file.path("/home", "steveped", "Documents", "Barry", "iClip", "alignments"), pattern=".sam",full.names=TRUE)
# out <- sam2bam(test, cl=cl, index=FALSE)
# stopCluster(cl)