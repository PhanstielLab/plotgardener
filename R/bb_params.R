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
#' p2 <- bb_params(fontsize = 10)
#'
#' ## Combine parameters and pass them to a BentoBox function
#' bb_plotGenes(params = c(p1, p2))
#'
#'
#' @export bb_params
#' @export c

## Define bb_params function skeleton (defined onLoad in zzz.R)
bb_params <- function(){}

## Define concatenate method for bb_params objects within default concatenate method to use when
## any bb_params object is found
"c" <- function(..., recursive = F){

  ## Check all classes of inputs to concatenate
  inputClasses <- unlist(lapply(list(...), class))
  ## If any are found to be `bb_params` objects, they will all be combined into one `bb_params` object
  if(any(inputClasses == "bb_params")){

    if (!all(inputClasses == "bb_params")){
      warning("Attempting to concatenate parameters not of class `bb_params` with `bb_params` objects. Coercing all parameters into a `bb_params` object.", call. = FALSE)

    }

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
                   paste(shQuote(badDupArgs), collapse = ",")), call. = FALSE)
    }

    ## Check that each allowed duplicate argument name has the same value
    for (a in allowed) {
      dups <- combArgs[names(combArgs) %in% a]
      if (length(unique(unlist(dups))) != 1) {
        stop(paste(sQuote(a), "must be the same when combining bb_params objects."), call. = FALSE)
      }
    }

    ## Return combined object
    return(structure(.Data = combArgs[!duplicated(names(combArgs))], class = "bb_params"))
  } else {
    ## Otherwise we will just call the primitive concatentate function
    .Primitive("c")(...)
  }

}
