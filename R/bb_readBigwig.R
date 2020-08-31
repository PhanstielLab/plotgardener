#' Reads a bigwig file and returns it as a data frame.
#'
#' @param filename bigwig filename
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param chrom chromosome as a string, e.g. "chr3"
#' @param chromstart first base number to consider
#' @param chromend last base number to consider
#' @param strand strand, i.e. '+', '-', or '*'
#'
#' @return Function will return a 6-column dataframe of bigwig information.
#' @export
#'
bb_readBigwig <- function(filename, params = NULL, chrom = NULL, chromstart = 1, chromend = .Machine$integer.max, strand = '*') {

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  bigwigDefaults <- c("chromstart", "chromend", "strand")

  ## Check which defaults are not overwritten and set to NULL
  if(missing(chromstart)) chromstart <- NULL
  if(missing(chromend)) chromend <- NULL
  if(missing(strand)) strand <- NULL

  ## Check if filename argument is missing (could be in object)
  if(!hasArg(filename)) filename <- NULL

  ## Compile all parameters into an internal object
  bb_bigwig <- structure(list(filename = filename, chrom = chrom, chromstart = chromstart, chromend = chromend, strand = strand), class = "bb_bigwig")

  bb_bigwig <- parseParams(bb_params = params, object_params = bb_bigwig)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_bigwig$chromstart)) bb_bigwig$chromstart <- 1
  if(is.null(bb_bigwig$chromend)) bb_bigwig$chromend <- .Machine$integer.max
  if(is.null(bb_bigwig$strand)) bb_bigwig$strand <- '*'

  # ======================================================================================================================================================================================
  # ERRORS
  # ======================================================================================================================================================================================

  if (is.null(bb_bigwig$filename)) stop("argument \"filename\" is missing, with no default.", call. = FALSE)

  if (bb_bigwig$chromend < bb_bigwig$chromstart - 1) {

    stop("End must be >= start - 1.", call. = FALSE)

  }

  # ======================================================================================================================================================================================
  # READ FILE
  # ======================================================================================================================================================================================

  if(!is.null(bb_bigwig$chrom)) {
    as.data.frame(rtracklayer::import.bw(bb_bigwig$filename, which=GenomicRanges::GRanges(paste0(bb_bigwig$chrom, ':', bb_bigwig$chromstart, '-', bb_bigwig$chromend, ':', bb_bigwig$strand))))
  } else {
    as.data.frame(rtracklayer::import.bw(bb_bigwig$filename))
  }
}
