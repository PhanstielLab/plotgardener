#' Reads a bigwig file and returns it as a data frame.
#'
#' @param filename bigwig filename
#' @param chrom chromosome as a string, e.g. "chr3"
#' @param chromstart first base number to consider
#' @param chromend last base number to consider
#' @param strand strand, i.e. '+', '-', or '*'
#' @export
#'
bb_readBigwig <- function(filename, chrom = NULL, chromstart = 1, chromend = .Machine$integer.max, strand = '*') {

  if (!hasArg(filename)) {

    stop("Filename must be specified.")

  }
  if (is.null(filename)) {

    stop("Filename cannot be NULL.")

  }
  if (chromend < chromstart - 1) {

    stop("End must be >= start - 1.")

  }


  if(!is.null(chrom)) {
    as.data.frame(import.bw(filename, which=GRanges(paste0(chrom, ':', chromstart, '-', chromend, ':', strand))))
  } else {
    as.data.frame(import.bw(filename))
  }
}
