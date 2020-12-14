#' bb_params: BentoBox parameters object
#'
#' Creates an object of class "bb_params" that can be used by BentoBox functions.
#' bb_params can be used to set a set of parameters to be shared across multiple
#' functions.
#'
#' bb_params generates arguments from exported BentoBox functions at loading time of the
#' package. Arguments defined in a bb_params object can be passed into the \code{params}
#' argument of BentoBox functions. bb_params arguments can be overridden from within
#' BentoBox functions.
#'
#' bb_params also provides an alternative region definition mechanism. Given a gene name
#' and genome assembly, bb_params returns the appropriate "chrom", "chromstart", and
#' "chromend" with a default buffer of \code{(gene length) / 2} added to the ends of the
#' gene coordinates. The buffer amount can be set manually with the \code{geneBuffer}
#' parameter. Buffer extending beyond the length of the chromosome will be trimmed.
#'
#' @param gene (optional) string naming a gene used to set the chromosome, chromstart, and
#'   chromend arguments.
#'
#' @param geneBuffer (optional) integer base-pairs to extend the start and end of a gene
#'   defined by argument \code{gene}.
#'
#' @param assembly string defining the genome build. Default value is \code{assembly = "hg19"}.
#'
#' @param ... This function will take any BentoBox function parameters and their values.
#'
#' @return Returns an object of class "bb_params" containing BentoBox function arguments.
#'
#' @examples
#' ## Define parameters
#' p1 <- bb_params(gene = "IL1B", assembly = "hg19")
#'
#' ## Optionally add more parameters
#' p2 <- bb_params(geneBuffer = 10000)
#'
#' ## Combine parameters and pass them to a BentoBox function
#' bb_plotGenes(params = c(p1, p2))
#'
#'
#' @export bb_params
#' @export c.bb_params

## Define bb_params function skeleton (defined onLoad in zzz.R)
bb_params <- function(){}

## Define concatenate method for bb_params objects
c.bb_params <- function(..., recursive = F) {

  ## Combine arguments into a single list
  combArgs <- unlist(list(...), recursive = F)

  ## Define allowed duplicates (i.e assembly="hg19")
  allowed <- c("assembly")

  ## Find duplicated argument names
  dupArgs <- combArgs[duplicated(names(combArgs))]

  ## Check for unallowed duplicate arguments
  if (any(!names(dupArgs) %in% allowed)) {
    badDupArgs <- names(dupArgs)[!names(dupArgs) %in% allowed]
    stop(sprintf("Parameter(s) %s are duplicated.",
                 paste(shQuote(badDupArgs), collapse = ",")))
  }

  ## Check that each allowed duplicate argument name has the same value
  for (a in allowed) {
    dups <- combArgs[names(combArgs) %in% a]
    if (length(unique(unlist(dups))) != 1) {
      stop(paste(sQuote(a), "must be the same when combining bb_params objects"))
    }
  }

  ## Return combined object
  return(structure(.Data = combArgs[!duplicated(names(combArgs))], class = "bb_params"))

}
