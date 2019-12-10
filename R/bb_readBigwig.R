#' Reads a bigwig file and returns it as a data frame.
#'
#' @param filename bigwig filename
#' @param chromosome chromosome as a string, e.g. "chr3"
#' @param start first base number to consider
#' @param end last base number to consider
#' @param strand strand, i.e. '+', '-', or '*'
#' @export
#'
bb_readBigwig <- function(filename, chromosome = NULL, start = 1, end = .Machine$integer.max, strand = '*') {

  if (!hasArg(filename)) {

    stop("Filename must be specified.")

  }
  if (is.null(filename)) {

    stop("Filename cannot be NULL.")

  }
  if (end < start - 1) {

    stop("End must be >= start - 1.")

  }


  if(!is.null(chromosome)) {
    as.data.frame(import.bw(filename, which=GRanges(paste0(chromosome, ':', start, '-', end, ':', strand))))
  } else {
    as.data.frame(import.bw(filename))
  }
}
