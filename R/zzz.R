.onLoad <- function(libname, pkgname) {

  ## Define bb_params object  -----------------------------

  ## Initialize documented bb_params function names
  x <- c("assembly", "gene", "geneBuffer")

  ## Set argument inputs for function definition
  allArgs1 <- paste(paste0(x, "=NULL"), collapse = ",")
  allArgs2 <- paste0(paste(x, "=", x), collapse = ",")

  ## Change specific argument defaults
  allArgs1 <- gsub('assembly=NULL', 'assembly="hg19"', allArgs1)
  allArgs1 <- paste0(allArgs1, ",...")

  ## Pass all arguments into function definition
  bb_params <- parse(text=c(sprintf("

  bb_params <- function(%s) {

    ## Construct object
    object <- structure(.Data = list(%s), class = 'bb_params')
    object[names(list(...))] <- list(...)

    ## Feature: setting region parameters by gene name & assembly ------------------------

    if (!is.null(gene)){

      ## Parse assembly
      assembly <- parse_bbAssembly(assembly = assembly)

      if (class(assembly$TxDb) == 'TxDb'){
        txdbChecks <- TRUE
      } else {
        txdbChecks <- check_loadedPackage(package = assembly$TxDb, message = paste(paste0('`', assembly$TxDb, '`'), 'not loaded. Please install and load to define genomic region based on a gene.'))
      }

      orgdbChecks <- check_loadedPackage(package = assembly$OrgDb, message = paste(paste0('`', assembly$OrgDb, '`'), 'not loaded. Please install and load to define genomic region based on a gene.'))

      if (txdbChecks == TRUE & orgdbChecks == TRUE){

        if (class(assembly$TxDb) == 'TxDb'){
          tx_db <- assembly$TxDb
        } else {
          tx_db <- eval(parse(text = assembly$TxDb))
        }

        org_db <- eval(parse(text = assembly$OrgDb))
        chromSizes <- GenomeInfoDb::seqlengths(tx_db)
        idCol <- assembly$gene.id.column
        displayCol <- assembly$display.column

        ## convert input gene name to geneID
        geneID <- tryCatch(expr = {
          suppressMessages(AnnotationDbi::select(org_db, keys = gene, columns = c(idCol), keytype = displayCol))

        }, error = function(e) stop(paste('Gene', shQuote(gene), 'does not exist in assembly.'), call. = FALSE))

        geneData <- suppressMessages(AnnotationDbi::select(tx_db, keys = geneID[[idCol]], columns = AnnotationDbi::columns(tx_db), keytype = 'GENEID'))


        ## Check that user has not supplied both gene and chrom, chromstart, or chromend
        chrom <- object$chrom
        chromstart <- object$chromstart
        chromend <- object$chromend
        if(any(!is.null(c(chrom, chromstart, chromend)))) {
          stop('Cannot use \\'gene\\' in combination with \\'chrom\\', \\'chromstart\\', or \\'chromend\\'', call. = FALSE)
        }

        ## Get info about gene region
        minGeneStart <- min(geneData$TXSTART)
        maxGeneEnd <- max(geneData$TXEND)

        ## Set default gene buffer (window = 2X gene length)
        ## Define buffer
        if (is.null(geneBuffer)) geneBuffer <- (maxGeneEnd - minGeneStart) / 2

        if (length(geneBuffer) == 1) geneBuffer <- rep(geneBuffer, 2)

        ## Assign values to bb_params object (with buffer)
        object$chrom      <- unique(geneData$TXCHROM)
        object$chromstart <- minGeneStart - geneBuffer[1]
        object$chromend   <- maxGeneEnd  + geneBuffer[2]
        object$geneBuffer <- geneBuffer

        ## Extract chromSizes length
        chrLength <- chromSizes[[object$chrom]]

        ## Check that starts and ends are within chromSizes
        if (object$chromstart < 1) {
          object$chromstart <- 1
          message('geneBuffer range is less than start. Start has been adjusted', call. = FALSE)
        }

        if (object$chromend > chrLength) {
          object$chromend   <- chrLength
          message('geneBuffer range is greater than end. End has been adjusted', call. = FALSE)
        }

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
