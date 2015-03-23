#' @title Index BAM files
#'
#' @description Indexes BAM files using \code{samtools index}
#' 
#' @details This is a wrapper to the function \code{samtools index} set to create the index files (.bai) to a supplied vector of bam files.
#' 
#' #' This function is designed for simple operations on files existing on a local HDD.
#' For more complex analyses or processes, please see the package \code{Rsamtools} on 
#' Bioconductor (\url{http://bioconductor.org/packages/release/bioc/html/Rsamtools.html}).
#' 
#' @param bamFiles A character vector with the full path of the BAM file(s)
#' @param cl A cluster as defined by the package \code{parallel}. If specified, the conversion will be run in parallel
#' @param exec The path to the samtools executable. Defaults to "/usr/bin/samtools". 
#' This may vary depending on where you have installed \code{samtools}
#' 
#' @return a list with the components \code{$bamFiles} & \code{$indexFiles}
#' 
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' 
#' @rdname samtoolsIndex
#' 
#' @import parallel
#' 
#' @export
samtoolsIndex <- function(bamFiles, cl=NULL, exec= "/usr/bin/samtools"){
  
  fEx <- file.exists(bamFiles)
  if (length(which(fEx)) == 0) stop("\nNo files were specified that exist")
  if (length(which(!fEx)) > 0 ) {
    warning("\nCould not find the following file(s), and they will be ignored:\n", bamFiles[!fEx] )
    bamFiles <- bamFiles[fEx]
  }
  
  areBam <- grep("[Bb][Aa][Mm]$", bamFiles)
  if (length(areBam) < length(bamFiles)) {
    warning("\nSome specified files appear to not be .bam files. These will be ignored\n")
    bamFiles <- bamFiles[areBam]
  }
  
  outFiles <- paste(bamFiles, ".bai", sep="")
  cmd <- paste(exec, "index", paste(bamFiles, outFiles))
  
  if (is.null(cl)) log <- sapply(cmd, FUN=system, intern=TRUE, USE.NAMES=FALSE)
  else log <- parSapply(cl, cmd, system, intern=TRUE, USE.NAMES=FALSE)
  fails <- which(unlist(log)==1)
  if (length(fails)>0) warning("\nThe indexing has failed on:\n", bamFiles[fails])
  
  return(list(bamFiles = bamFiles, indexFiles = outFiles))
  
}

# cl <- makeCluster(8, type = "SOCK")
# test <- list.files(file.path("/home", "steveped", "Documents", "Barry", "iClip", "alignments"), pattern=".bam",full.names=TRUE)
# out <- samtoolsIndex(test, cl=cl)
# stopCluster(cl)