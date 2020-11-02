#' Creates a bb_params object that can be used by BentoBox functions
#'
#' @param gene Optional parameter to extract chromosome, chromstart, and chromend information for a gene
#' @param assembly If specifying gene, default genome assembly as a string or a bb_assembly object
#' @param ... This function will take any BentoBox function parameters and their values
#'
#' @return All function arguments will be returned as a "bb_params" object
#' @export

bb_params <- function(gene = NULL, assembly = "hg19",...){

  ## Masterlist of all possible BentoBox parameters to check that it falls into that and throws warning if not

  object <- structure(list(...), class = "bb_params")
  if (!is.null(gene)){
    object$gene <- gene
    ### CHANGE CODE HERE TO GET ASSOCIATED TXDB/ORG FOR GENE
    if (assembly == "hg19"){
      genes <- bb_hg19gtf
    }
    geneRegion <- genes[which(genes$Gene == bb_params$gene),]
    object$chrom <- geneRegion$Chromosome
    object$chromstart <- geneRegion$Start
    object$chromend <- geneRegion$Stop

  }


  return(object)

}
