#' @title Read a single log file as written by hisat2
#'
#' @description Takes a log file as created by hisat2, and loads as a data_frame with a simgle row
#' 
#' @details Loads in the log file and extracts the key values to form a summary.
#' Multiple log Files can be loaded using \code{lapply} followed by \code{bind_rows}
#' 
#' @param file The path to a single file
#' 
#' @return A \code{\link{data_frame}} or \code{\link{tbl_df}}
#' 
#' @author Steve Pederson <stephen.pederson@@adelaide.edu.au>
#'
#' @import dplyr
#' @export
readHisat2Log <- function(file){
  if (length(file) > 1) {
    warning("Multiple files specified. Only the first will be read. Please use lapply for multiple files")
    file <- file[1]
  }
  stopifnot(file.exists(file))
  ln <- readLines(file)
  
  dplyr::data_frame(
    File = basename(file),
    TotalReads = as.integer(gsub("([0-9]*) reads; of these:", "\\1", ln[1])),
    PairedReads = as.integer(gsub("([0-9]*) \\(.+\\) were paired; of these:", "\\1", ln[2])),
    UniqueInPairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned concordantly exactly 1 time", "\\1", ln[4])),
    MultipleInPairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned concordantly >1 times", "\\1", ln[5])),
    UniqueDiscordantPairs = as.integer(gsub("([0-9]*) \\(.+\\) aligned discordantly 1 time", "\\1", ln[8])),
    UniqueUnpaired = as.integer(gsub("([0-9]*) \\(.+\\) aligned exactly 1 time", "\\1", ln[13])),
    MultipleUnpaired = as.integer(gsub("([0-9]*) \\(.+\\) aligned >1 times", "\\1", ln[14])),
    NotAligned = as.integer(gsub("([0-9]*) \\(.+\\) aligned 0 times", "\\1", ln[12])),
    AlignmentRate = 1 - NotAligned / (TotalReads + PairedReads)
  )
}