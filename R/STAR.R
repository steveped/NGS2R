#' @title A wrapper for the RNA-Seq Aligner STAR
#'
#' @description Builds a STAR index, or align to the genome using STAR.
#' 
#' @details \code{buildSTARIndex} will build the index for aligning using the STAR RNA-seq aligner. 
#' This tool is built for multi-threading and will run much faster if this capbility is taken advantage of.
#' 
#' \code{alignSTAR} is the wrapper for aligning a single fastq file, or a single set of paired end reads.
#' For multiple samples, please use the additional wrappers \code{lapply} or \code{snow::clusterApply}.
#' However, if aligning in parallel, please set the parameter \code{threads=1} as the function will automatically use half of the available cores, as STAR itself is written wil mutli-threading capabilities.
#' 
#' The commands for correctly utilising the shared memory capabilities of STAR have been problematic to implement, 
#' so if you specifically require these, please use STAR on the command line.
#' 
#' @param genome The reference genome as a single fasta file.
#' @param indexPath The directory to write the index files to when bilding, or to look for the index files when aligning
#' @param threads The number of threads to utilise. Defults to half the available cores.
#' @param gtfFile The annotation file in \code{.gtf} format
#' @param overhang The length of the overhang permitted. This should be set as max(ReadLength) - 1. Will default to 84, but this may not be suitable so please check this parameter
#' @param inFiles The input file(s). Should be a single fastq file, or paired end reads from a single sample.
#' @param outDir The directory to write the alignments to.
#' @param outType Specify whether to write the output files as SAM or BAM format. Defaults to BAM
#' @param sorted Specify whether to sort the alignments (\code{TRUE}) in the output files. Defaults to true.
#' @param prefix Any prefix to add to the sam/bam files. Will default to the prefix on the fastq file(s)
#' @param args Any additional arguments to pass to STAR. 
#' Please note these will not be checked for errors, so please check the STAR manual carefully \url{https://github.com/alexdobin/STAR}.
#' @param exec The path to the STAR executable
#' 
#' @return All functions execute the tool on the command line, then return a \code{list} with key information
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' 
#' @import parallel
#' 
#' @rdname STAR
#' @export

buildSTARIndex <- function(genome, indexPath, threads, gtfFile, overhang=84, exec="/usr/local/bin/STAR"){
  
  if (!file.exists(genome)) stop("\nCannot find specified genome files:\n", genome)
  if (!gsub(".+\\.([A-Za-z0-9]+$)","\\1", basename(genome)) %in% c("fa", "mfa", "fasta")) stop("")
  
  
  if (!file.exists(indexPath)) dir.create(indexPath, recursive=TRUE)
  
  nCores <- detectCores()
  if (missing(threads)) threads <- nCores/2
  threads <- as.integer(min(threads, nCores)) # Ensure you can't have more threads than cores
  
  if (!missing(gtfFile)) {
    if (!file.exists(gtfFile)) stop("\nCannot find specified gtf file:\n", gtf)
    if (gsub(".+\\.([A-Za-z0-9]+$)","\\1", basename(gtfFile)) != "gtf") stop("\nExpected suffix .gtf:\n", gtfFile)
    gtf <- paste("--sjdbGTFfile", gtfFile)
  }
  else {
    gtf <- c()
  }
  
  cmd <- paste(exec,
               "--runMode genomeGenerate", 
               "--runThreadN", threads,
               "--genomeDir", indexPath,
               "--genomeFastaFiles", genome,
               gtf,
               "--sjdbOverhang", overhang)
  
  st <- Sys.time()
  cat("Commencing index build at", date(), "\n")
  
  vers <- system(paste(exec, "--version"), intern=TRUE)
  cat("Building index using:", vers[1], "\n")

  system(cmd, intern=TRUE)
  ft <- Sys.time()
  cat("Index build complete at", date(), "\n")
 
  return(list(indexPath = indexPath, version = vers, command = cmd, time = data.frame(start=st, finish=ft)))
  
}

#' @rdname STAR
#' @export
alignSTAR <- function(inFiles, indexPath, outDir, outType="BAM", sorted=TRUE, threads, prefix, args = c(), exec = "/usr/local/bin/STAR"){
  
  fe <- which(!file.exists(inFiles))
  if (length(fe) > 0) stop("\nCannot find file(s):\n", inFiles[fe])
  
  if (!file.exists(indexPath)) 
    stop("\nCannot find index:\n", indexPath)
  indexFiles <- list.files(indexPath)
  starFiles <- c(paste(c("chrLength", "chrName", "chrNameLength", "chrStart", "genomeParameters", "sjdbInfo"), 
                       "txt", sep="."),
                 paste(c("exonGeTrInfo", "exonInfo", "geneInfo", "sjdbList.fromGTF.out", "sjdbList.out", "transcriptInfo"),
                       "tab", sep="."),
                 c("Genome", "SA", "SAindex"))
  if (length(which(!starFiles %in% indexFiles)) > 0) stop("\nIt looks like some of the required files for the index are missing\n")
  cat("Genome index files found in", indexPath, "\n")
  
  if (!file.exists(outDir)) dir.create(outDir, recursive=TRUE)
  
  if (length(outType)>1) warning("More than one value was set for outType. Only the first will be used and other will be ignored\n", immediate.=TRUE)
  outType <- toupper(outType[1])
  if (!outType %in% c("SAM", "BAM")) stop("\nOutput must be specified as either BAM or SAM")
  outCmd <- ifelse(sorted, paste(outType, "SortedByCoordinate"), paste(outType, "Unsorted"))
  
  nCores <- detectCores()
  if(missing(threads) || !is.numeric(threads)) threads <- nCores/2
  threads <- as.integer(min(threads, nCores))
  
  
  if (missing(prefix)) prefix <- gsub("R[12]", "", 
                                      gsub("\\.(fq|fastq|fq\\.gz|fastq\\.gz)", "", basename(inFiles)))[1]
  outPrefix <- file.path(outDir, prefix)
  
  ## Now pull the command together & run...
  cmd <- paste(exec, "--genomeDir", indexPath,
               "--runThreadN", threads,
               "--readFilesIn", paste(inFiles, collapse=" "),
               "--outFileNamePrefix", outPrefix, 
               "--outSAMtype", outCmd,
               args)
  
  vers <- system("STAR --version", intern=TRUE)
  
  cat("Using STAR version", vers, "\n")
  cat("Executing system command:\n", cmd, "\n")
  system(cmd)
  
  pat <- paste(prefix, ".+", tolower(outType),"$", sep="")
  alnFiles <- list.files(outDir, full.names=TRUE)
  return(list(alnFiles = grep(pat, alnFiles, value=TRUE), 
              logFiles = grep("Log", alnFiles, value=TRUE),
              tabFiles = grep("tab", alnFiles, value=TRUE),
              outDir = outDir,
              version = vers,
              cmd=cmd))
  
}
