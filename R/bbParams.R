#' bbParams: BentoBox parameters object
#'
#' Creates an object of class "bbParams" that can be used by
#' BentoBox functions. bbParams can be used to set a set of parameters
#' to be shared across multiple functions.
#'
#' bbParams generates arguments from exported BentoBox functions at
#' loading time of the package. Arguments defined in a bbParams object
#' can be passed into the \code{params} argument of BentoBox functions.
#' bbParams arguments can be overridden from within BentoBox functions.
#'
#' bbParams also provides an alternative region definition mechanism.
#' Given a gene name and genome assembly, bbParams returns the appropriate
#' "chrom", "chromstart", and "chromend" with a default buffer of
#' \code{(gene length) / 2} added to the ends of the gene coordinates.
#' The buffer amount can be set manually with the \code{geneBuffer}
#' parameter. Buffer extending beyond the length of the chromosome
#' will be trimmed.
#'
#' @param gene (optional) String naming a gene used to set the
#' chrom, chromstart, and chromend arguments.
#'
#' @param geneBuffer (optional) Integer base-pairs to extend the
#' start and end of a gene defined by argument \code{gene}.
#' Can be one integer or a vector of length 2, where the first integer
#' will extend the start of the gene and the second integer
#' will extend the end of the gene.
#'
#' @param assembly String defining the genome build.
#' Default value is \code{assembly = "hg38"}.
#' 
#'
#' @param ... This function will take any BentoBox function
#' parameters and their values:
#' \itemize{
#' \item{\code{alpha}}
#' \item{\code{altchrom}}
#' \item{\code{altchromend}}
#' \item{\code{altchromstart}}
#' \item{\code{archHeight}}
#' \item{\code{arrow}}
#' \item{\code{at}}
#' \item{\code{axis}}
#' \item{\code{axisLine}}
#' \item{\code{baseline}}
#' \item{\code{baseline.color}}
#' \item{\code{baseline.lwd}}
#' \item{\code{bg}}
#' \item{\code{binCap}}
#' \item{\code{binSize}}
#' \item{\code{border}}
#' \item{\code{boxHeight}}
#' \item{\code{boxWidth}}
#' \item{\code{breaks}}
#' \item{\code{BSgenome}}
#' \item{\code{cex}}
#' \item{\code{check.overlap}}
#' \item{\code{chrom}}
#' \item{\code{chromend}}
#' \item{\code{chromstart}}
#' \item{\code{clip}}
#' \item{\code{collapse}}
#' \item{\code{colorbyStrand}}
#' \item{\code{colorTrans}}
#' \item{\code{column}}
#' \item{\code{commas}}
#' \item{\code{curvature}}
#' \item{\code{data}}
#' \item{\code{default.units}}
#' \item{\code{digits}}
#' \item{\code{display.column}}
#' \item{\code{draw}}
#' \item{\code{extend}}
#' \item{\code{file}}
#' \item{\code{fill}}
#' \item{\code{flip}}
#' \item{\code{fontcolor}}
#' \item{\code{fontsize}}
#' \item{\code{geneBackground}}
#' \item{\code{geneHighlights}}
#' \item{\code{gene.id.column}}
#' \item{\code{geneOrder}}
#' \item{\code{Genome}}
#' \item{\code{half}}
#' \item{\code{height}}
#' \item{\code{id}}
#' \item{\code{id.lengths}}
#' \item{\code{image}}
#' \item{\code{interpolate}}
#' \item{\code{just}}
#' \item{\code{label}}
#' \item{\code{labels}}
#' \item{\code{leadSNP}}
#' \item{\code{legend}}
#' \item{\code{length}}
#' \item{\code{limitLabel}}
#' \item{\code{linecolor}}
#' \item{\code{lineend}}
#' \item{\code{linejoin}}
#' \item{\code{lty}}
#' \item{\code{lwd}}
#' \item{\code{main}}
#' \item{\code{margin}}
#' \item{\code{matrix}}
#' \item{\code{negData}}
#' \item{\code{norm}}
#' \item{\code{OrgDb}}
#' \item{\code{orientation}}
#' \item{\code{palette}}
#' \item{\code{pch}}
#' \item{\code{plot}}
#' \item{\code{quiet}}
#' \item{\code{r}}
#' \item{\code{range}}
#' \item{\code{resolution}}
#' \item{\code{res_scale}}
#' \item{\code{rot}}
#' \item{\code{scale}}
#' \item{\code{scientific}}
#' \item{\code{scipen}}
#' \item{\code{sequence}}
#' \item{\code{shift}} 
#' \item{\code{showBands}}
#' \item{\code{showGuides}}
#' \item{\code{sigCol}}
#' \item{\code{sigLine}}
#' \item{\code{sigVal}}
#' \item{\code{spaceHeight}}
#' \item{\code{spaceWidth}}
#' \item{\code{strand}}
#' \item{\code{strandLabels}}
#' \item{\code{strandSplit}}
#' \item{\code{stroke}}
#' \item{\code{style}}
#' \item{\code{tcl}}
#' \item{\code{ticks}}
#' \item{\code{title}}
#' \item{\code{TxDb}}
#' \item{\code{type}}
#' \item{\code{width}}
#' \item{\code{x}}
#' \item{\code{xgrid}}
#' \item{\code{x0}}
#' \item{\code{x1}}
#' \item{\code{y}}
#' \item{\code{ygrid}}
#' \item{\code{ymax}}
#' \item{\code{y0}}
#' \item{\code{y1}}
#' \item{\code{zrange}}
#' }
#'
#' @return Returns an object of class \code{bbParams}
#' containing BentoBox function arguments.
#'
#' @examples
#' ## Load hg19 genomic annotation packages
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Define parameters
#' p1 <- bbParams(gene = "IL1B", assembly = "hg19")
#'
#' ## Optionally add more parameters
#' p2 <- bbParams(fontsize = 10, assembly = "hg19")
#'
#' ## Combine parameters and pass them to a BentoBox function
#' bbPlotGenes(params = c(p1, p2))
#' @export bbParams
#' @export c

## Define bbParams function skeleton (defined onLoad in zzz.R)
bbParams <- function() {}

## Define concatenate method for bbParams objects within default concatenate
## method to use when
## any bbParams object is found

#' Combine multiple bbParams objects into a vector
#'
#' @param ... \link[BentoBox]{bbParams} objects to be concatenated.
#' @param recursive logical. If \code{recursive = TRUE}, the function
#' recursively descends through lists
#' (and pairlists) combining all their elements into a vector.
#'
#' @return \code{NULL} or an expression or a vector of an appropriate mode.
#' (With no arguments the value is \code{NULL}.)
#'
#' @examples
#' ## Define parameters
#' p1 <- bbParams(chrom = "chr1", assembly = "hg19")
#'
#' ## Define another set of parameters
#' p2 <- bbParams(fontsize = 10, assembly = "hg19")
#'
#' ## Combine parameters into one `bbParams` object
#' pTotal <- c(p1, p2)
"c" <- function(..., recursive = FALSE) {

    ## Check all classes of inputs to concatenate
    inputClasses <- unlist(lapply(list(...), class))
    ## If any are found to be `bbParams` objects, they will all be combined
    ## into one `bbParams` object
    if (any(inputClasses == "bbParams")) {
        if (!all(inputClasses == "bbParams")) {
            warning("Attempting to concatenate parameters not of ",
                    "class `bbParams` with `bbParams` objects. ",
                    "Coercing all parameters into a `bbParams` object.",
                call. = FALSE
            )
        }

        ## Combine arguments into a single list
        combArgs <- unlist(list(...), recursive = FALSE)
        ## Define allowed duplicates (i.e assembly="hg19")
        allowed <- c("assembly")

        ## Find duplicated argument names
        dupArgs <- combArgs[duplicated(names(combArgs))]

        ## Check for duplicate arguments that aren't allowed
        if (any(!names(dupArgs) %in% allowed)) {
            badDupArgs <- names(dupArgs)[!names(dupArgs) %in% allowed]
            message <- sprintf(
                "Parameter(s) %s are duplicated.",
                paste(shQuote(badDupArgs), collapse = ",")
            )
            stop(message, call. = FALSE)
        }

        ## Check that each allowed duplicate argument name has the same value
        for (a in allowed) {
            dups <- combArgs[names(combArgs) %in% a]
            if (length(unique(unlist(dups))) != 1) {
                stop(sQuote(a), "must be the same when combining ",
                    "bbParams objects.", call. = FALSE)
            }
        }

        ## Return combined object
        return(structure(
            .Data = combArgs[!duplicated(names(combArgs))],
            class = "bbParams"
        ))
    } else {
        ## Otherwise we will just call the primitive concatenate function
        .Primitive("c")(...)
    }
}
