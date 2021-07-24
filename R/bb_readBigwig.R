#' Read a bigWig file and return it as a data frame
#' 
#' @usage bb_readBigwig(
#'     file,
#'     chrom = NULL,
#'     chromstart = 1,
#'     chromend = .Machine$integer.max,
#'     strand = "*",
#'     params = NULL
#' )
#'
#' @param file A character value specifying the path to the bigwig file.
#' @param chrom Chromosome of data as a string, if data for a specific
#' chromosome is desired.
#' @param chromstart Integer start position on chromosome.
#' @param chromend Integer end position on chromosome.
#' @param strand A character value specifying strand.
#' Default value is \code{strand = "*"}. Options are:
#' \itemize{
#' \item{\code{"+"}: }{Plus strand.}
#' \item{\code{"-"}: }{Minus strand.}
#' \item{\code{"*"}: }{Plus and minus strands.}
#' }
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#'
#' @return Returns a 6-column dataframe of bigwig information.
#' 
#' @examples 
#' if (.Platform$OS.type != "windows"){
#'     bwFile <- system.file("extdata/test.bw", package="BentoBoxData")
#' 
#'     ## Read in entire file
#'     bwData <- bb_readBigwig(file = bwFile)
#' 
#'     ## Read in specified region
#'     bwRegion <- bb_readBigwig(file = bwFile,
#'                             chrom = "chr2",
#'                             chromstart = 1,
#'                             chromend = 1500)
#' }
#' 
#' @details This function does not work on Windows.
#' 
#' @seealso \link[rtracklayer]{import.bw}
#' 
#' @export
bb_readBigwig <- function(file, chrom = NULL, chromstart = 1,
                        chromend = .Machine$integer.max,
                        strand = "*", params = NULL) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================
        
    bb_bigwig <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_bigwig"
    )

    # =========================================================================
    # ERRORS
    # =========================================================================

    if (is.null(bb_bigwig$file)) {
        stop("argument \"file\" is missing, with no default.",
            call. = FALSE
        )
    }

    if (bb_bigwig$chromend < bb_bigwig$chromstart - 1) {
        stop("End must be >= start - 1.", call. = FALSE)
    }

    # =========================================================================
    # READ FILE
    # =========================================================================

    if (!is.null(bb_bigwig$chrom)) {
        as.data.frame(rtracklayer::import.bw(bb_bigwig$file,
            which = GenomicRanges::GRanges(paste0(
                bb_bigwig$chrom,
                ":",
                bb_bigwig$chromstart,
                "-",
                bb_bigwig$chromend,
                ":",
                bb_bigwig$strand
            ))
        ))
    } else {
        as.data.frame(rtracklayer::import.bw(bb_bigwig$file))
    }
}
