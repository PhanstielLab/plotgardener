.onLoad <- function(libname, pkgname) {

  ## Define bb_params object with all exported function args -----------------------------

  ## Extract exported function names
  fname <- paste0(system.file(lib.loc = libname, package = pkgname),"/NAMESPACE")
  fxns <- gsub("export\\((.*)\\)", "\\1", grep("export", readLines(fname), value = T))

  ## Get arguments from all functions
  x <- unlist(lapply(fxns, function(x) names(formals(x))))

  ## Filter out arguments
  x <- x[!x %in% c("", "...", "params", "recursive")]

  ## Add in bb_params-specific arguments (do this here for correct sorting, default=NULL)
  x <- unique(c("gene", "geneBuffer", x))

  ## Alphebetize
  x <- sort(x, decreasing = F)

  ## Set argument inputs for function definition
  allArgs1 <- paste(paste0(x, "=NULL"), collapse = ",")
  allArgs2 <- paste0(paste(x, "=", x), collapse = ",")

  ## Change specific argument defaults (add lines as desired ...)
  allArgs1 <- gsub('assembly=NULL', 'assembly="hg19"', allArgs1)

  ## Pass all arguments into function definition
  bb_params <- parse(text=c(sprintf("

  bb_params <- function(%s) {

    ## Construct object
    object <- structure(.Data = list(%s), class = 'bb_params')

    ## Feature: setting region parameters by gene name & assembly ------------------------

    if (!is.null(gene)){

      ## CHANGE CODE HERE TO GET ASSOCIATED TXDB/ORG FOR GENE
      if (assembly == 'hg19')
      {
        genes <- bb_hg19gtf
        chromSizes <- bb_hg19
      }
      else if (assembly == 'hg38')
      {
        genes <- bb_hg38gtf
        chromSizes <- bb_hg38
      }
      else
      {
        stop(paste('Assembly', shQuote(assembly), 'is not defined.'))
      }

      ## Check that gene is in gtf file
      if(!gene %%in%% genes$Gene) {
        stop(paste('Gene', shQuote(gene),
                   'does not exist in assembly', shQuote(assembly)))
      }

      ## Check that user has not supplied both gene and chrom, chromstart, or chromend
      if(any(!is.null(c(chrom, chromstart, chromend)))) {
        stop('Cannont use \\'gene\\' in combination with \\'chrom\\', \\'chromstart\\', or \\'chromend\\'')
      }

      ## Subset for gene region
      geneRegion <- genes[genes$Gene == object$gene,]


      ## Set default gene buffer (window = 2X gene length)

      ## Define buffer
      if (is.null(geneBuffer)) geneBuffer <- (geneRegion$Stop - geneRegion$Start) / 2

      ## Assign values to bb_params object (with buffer)
      object$chrom      <- geneRegion$Chromosome
      object$chromstart <- geneRegion$Start - geneBuffer
      object$chromend   <- geneRegion$Stop  + geneBuffer
      object$geneBuffer <- geneBuffer

      ## Extract chromSizes length
      chrLength <- chromSizes$length[chromSizes$chrom == object$chrom]

      ## Check that starts and ends are within chromSizes
      if (object$chromstart < 1) {
        object$chromstart <- 1
        message('geneBuffer range is less than start. Start has been adjusted')
      }
      if (object$chromend > chrLength) {
        object$chromend   <- chrLength
        message('geneBuffer range is greater than end. End has been adjusted')
      }

    }

    ## Filter out null values for pretty printing
    object <- structure(Filter(Negate(is.null), object), class = 'bb_params')

    return(object)

  }

  ", allArgs1, allArgs2)))

  ## Assign function in BentoBox environment
  assign("bb_params", eval(bb_params), rlang::ns_env(pkgname))

}
