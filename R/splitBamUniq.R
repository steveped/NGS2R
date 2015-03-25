#' @title Split a BAM file into unique alignments & multiple alignments
#'
#' @description Creates two new BAM files based on a single BAM file
#'
#' @details This is a wrapper to the function \code{samtools view} and \code{grep}.
#' 
#' Files will be split and the names of the newly created bam files will be supplied in the output in the component \code{$outFiles}.
#' Currently the conversion takes place via SAM files, which will be converted back to BAM files.
#' The SAM files will then be deleted.
#' 
#' This function is designed for simple operations on files existing on a local HDD.
#' For more complex analyses or processes, please see the package \code{Rsamtools} on 
#' Bioconductor (\url{http://bioconductor.org/packages/release/bioc/html/Rsamtools.html})
#'
#' @param inFile The BAM file which is to be split
#' @param outPath The directory to write to. Defaults to the same directory as \code{inFile}. Will be created if it doesn't already exist.
#' @param exec The path to the samtooles executable 
#'
#' @return A list with the components \code{inFile} & \code{outFiles}
#'
#' @author Steve Pederson
#' 
#' @rdname splitBamUniq
#' 
#' @export
splitBamUniq <- function(inFile, outPath, exec = "/usr/bin/samtools"){
  
  if (length(inFile) > 1) {
    warning("More than one file supplied. Only the first will be split & all other will be ignored.\n", immediate.=TRUE)
    inFile <- inFile[1]
  }
  if (!file.exists(inFile)) stop("Cannot find the requested BAM file. Please check you have the correct path and file name.\n")
  isBam <- as.logical(regexpr("[B|b][A|a][M|m]$", inFile)[1] + 1)
  if (!isBam) stop("The supplied BAM file must have the suffix .bam\n")
  
  if (missing(outPath)) {
    outPath <- dirname(inFile)
  }
  else {
    if (!file.exists(outPath))  dir.create(outPath, recursive=TRUE)
  }
  
  if (!file.exists(exec)) stop("Cannot find the executable. Please check the location inthe terminal using the command whereis\n")
    
  outSam <- list(unique = gsub(".[B|b][A|a][M|m]$", ".unique.sam", basename(inFile)),
                 multi = gsub(".[B|b][A|a][M|m]$", ".multi.sam", basename(inFile)))
  outSam <- lapply(outSam, function(x){file.path(outPath, x)})
    
  cloneHeader(inFile, outSam$unique, exec=exec)
  cloneHeader(inFile, outSam$multi, exec=exec)
  fEx <- sapply(outSam, file.exists)
  if (any(!fEx)) stop("There seems to have been a problem creating the initial SAM files. Do you have permission to write to that directory?\n")
  
  addUniqAlignments(inFile, outSam$unique, exec)
  addMultiAlignments(inFile, outSam$multi, exec)
  
  out <- sam2bam(unlist(outSam), sort=FALSE, exec = exec)
  fEx <- file.exists(out$unsorted)
  if (any(!fEx)) stop("There seems to have been a problem converting the SAM files to BAM files.\n")
  
  file.remove(out$inFiles)
  return(list(inFile = inFile, outFiles = out$unsorted))
  
}

cloneHeader <- function(inFile, outFile, exec = "/usr/bin/samtools"){
  type <- grep(tolower(gsub(".+\\.(.{3})$", "\\1", inFile)), c("sam", "bam"), value=TRUE)
  if (type=="bam") cmd <- paste(exec, "view -H", inFile, ">", outFile)
  else cmd <- paste(exec, "view -H -S", inFile, ">", outFile)
  system(cmd, intern=TRUE)
  return(outFile)
}

addUniqAlignments <- function(inFile, outFile, exec = "/usr/bin/samtools"){
  cmd <- paste(exec, "view", inFile, "| grep 'NH:i:1\\s' >>", outFile)
  system(cmd)
  return(outFile)
}

addMultiAlignments <- function(inFile, outFile, exec = "/usr/bin/samtools"){
  cmd <- paste(exec, "view", inFile, "| grep -v 'NH:i:1\\s' >>", outFile)
  system(cmd)
  return(outFile)
}