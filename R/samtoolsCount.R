#' @title Count the requested number of reads in a given BAM/SAM File
#'
#' @description Uses the requested filters to provide the count of reads which match within the specified BAM or SAM file.
#' 
#' @details This is a wrapper to the function \code{samtools view} with the option \code{-c} set.
#' Any combination of flags or quality filters can be used to extract the requested counts.
#' 
#' This function is designed for simple operations on files existing on a local HDD.
#' For more complex analyses or processes, please see the package \code{Rsamtools} on 
#' Bioconductor (\url{http://bioconductor.org/packages/release/bioc/html/Rsamtools.html}).
#' 
#' @param file A single BAM/SAM file to query
#' @param include Sets the flags for inclusion. This is the command-line parameter "-f"
#' @param exclude Sets the flags for exclusion. This is the command-line parameter "-F"
#' @param minQ Sets the minimum mapping quality to include
#' @param unique Provide the number of unique reads which match the specified criteria. Defaults to FALSE
#' @param exec The path to the samtools executable. Defaults to "/usr/bin/samtools". 
#' This may vary depending on where you have installed \code{samtools}
#' 
#' @return an integer corresponding to the count of reads matching the specified filter
#' 
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#' 
#' @rdname samtoolsCount
#' 
#' @export
samtoolsCount <- function(file, include=NULL, exclude=NULL, minQ=NULL, unique=FALSE, exec = "/usr/bin/samtools"){
  
  if (length(file) > 1) warning("Too many files. Only the firsty file will be processed.\n")
  file <- file[1]
  if (!file.exists(file)) stop("\nCannot find file:\n", file)
  type <- tolower(gsub("^.+\\.([BSbs][Aa][Mm]$)","\\1", basename(file)))
  if(type=="bam") {
    S <- c()
  }
  else{
    if(type=="sam") S <- "-S"
    else stop("\nUnknown file type. Only .bam or .sam files can be requested\n")
  }
  
  if (!is.null(include)) {
    if (length(include) > 1 ) warning("Only single integers can be accepted for inclusion. Only the first value will be used")
    if (!is.numeric(include)) stop("\nFlags for inclusion can only be specified numerically:", include)
    include <- paste("-f", as.integer(include))
  }
  if (!is.null(exclude)) {
    if (length(exclude) > 1 ) warning("Only single integers can be accepted for exclusion. Only the first value will be used")
    if (!is.numeric(exclude)) stop("\nFlags for exclusion can only be specified numerically:", exclude)
    exclude <- paste("-F", as.integer(exclude))
  }
  if (!is.null(minQ)) {
    if (length(minQ) > 1 ) warning("Only single integers can be accepted for quality filtering. Only the first value will be used")
    if (!is.numeric(minQ)) stop("\nQuality filters can only be specified numerically:", minQ)
    minQ <- paste("-q", as.integer(minQ))
  }

  if (!is.logical(unique)) stop("\nThe parameter unique must be specified as a logical value\n")
  unique <- unique[1]
  
  if (!unique) cmd <- paste(exec, "view -c", S, include, exclude, minQ, file)
  else cmd <- paste(exec, "view", S, include, exclude, minQ, file, "| cut -f1 | sort | uniq | wc -l")
  cmd <- gsub(" +", " ", cmd)
  return(as.integer(system(cmd, intern=TRUE, ignore.stderr=TRUE)))

}
