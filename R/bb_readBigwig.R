#' Read a bigwig file and return it as a data frame
#'
#' @usage
#' bb_readBigwig(bigwigFile)
#' bb_readBigwig(bigwigFile, chrom, chromstart, chromend)
#'
#' @param bigwigFile A character value specifying the path to the bigwig file.
#' @param chrom Chromosome of data as a string, if data for a specific chromosome is desired.
#' @param chromstart Integer start position on chromosome.
#' @param chromend Integer end position on chromosome.
#' @param strand A character value specifying strand. Default value is \code{strand = "*"}.Options are:
#' \itemize{
#' \item{\code{"+"}: }{Plus strand.}
#' \item{\code{"-"}: }{Minus strand.}
#' \item{\code{"*"}: }{Plus and minus strands.}
#' }
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#'
#' @return Returns a 6-column dataframe of bigwig information.
#'
#' @export
bb_readBigwig <- function(bigwigFile, chrom = NULL, chromstart = 1, chromend = .Machine$integer.max, strand = '*', params = NULL) {

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  bigwigDefaults <- c("chromstart", "chromend", "strand")

  ## Check which defaults are not overwritten and set to NULL
  if(missing(chromstart)) chromstart <- NULL
  if(missing(chromend)) chromend <- NULL
  if(missing(strand)) strand <- NULL

  ## Check if filename argument is missing (could be in object)
  if(!hasArg(bigwigFile)) bigwigFile <- NULL

  ## Compile all parameters into an internal object
  bb_bigwig <- structure(list(bigwigFile = bigwigFile, chrom = chrom, chromstart = chromstart, chromend = chromend, strand = strand), class = "bb_bigwig")

  bb_bigwig <- parseParams(bb_params = params, object_params = bb_bigwig)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_bigwig$chromstart)) bb_bigwig$chromstart <- 1
  if(is.null(bb_bigwig$chromend)) bb_bigwig$chromend <- .Machine$integer.max
  if(is.null(bb_bigwig$strand)) bb_bigwig$strand <- '*'

  # ======================================================================================================================================================================================
  # ERRORS
  # ======================================================================================================================================================================================

  if (is.null(bb_bigwig$bigwigFile)) stop("argument \"bigwigFile\" is missing, with no default.", call. = FALSE)

  if (bb_bigwig$chromend < bb_bigwig$chromstart - 1) {

    stop("End must be >= start - 1.", call. = FALSE)

  }

  # ======================================================================================================================================================================================
  # READ FILE
  # ======================================================================================================================================================================================

  if(!is.null(bb_bigwig$chrom)) {
    as.data.frame(rtracklayer::import.bw(bb_bigwig$bigwigFile, which=GenomicRanges::GRanges(paste0(bb_bigwig$chrom, ':', bb_bigwig$chromstart, '-', bb_bigwig$chromend, ':', bb_bigwig$strand))))
  } else {
    as.data.frame(rtracklayer::import.bw(bb_bigwig$bigwigFile))
  }
}
