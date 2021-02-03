#' Read a .hic file and return Hi-C data as a dataframe
#'
#' @param file A character value specifying the path to the .hic file.
#' @param chrom Chromosome of data, as a string.
#' @param chromstart Integer start position on chromosome.
#' @param chromend Integer end position on chromosome.
#' @param altchrom Alternate chromosome for interchromosomal data, as a string.
#' @param altchromstart Alternate chromosome integer start position for interchromosomal data.
#' @param altchromend Alternate chromosome integer end position for interchromosomal data.
#' @param assembly Default genome assembly as a string or a \link[BentoBox]{bb_assembly} object. Default value is \code{assembly = "hg19"}.
#' @param resolution A numeric specifying the width of each pixel. "auto" will attempt to choose a resolution in basepairs based on the size of the region.
#' @param res_scale A character value specifying the resolution scale. Default value is \code{res_scale = "BP"}. Options are:
#' \itemize{
#' \item{\code{"BP"}: }{Base pairs.}
#' \item{\code{"FRAG"}: }{Fragments.}
#' }
#' @param zrange A numeric vector of length 2 specifying the range of interaction scores, where extreme values will be set to the max or min.
#' @param norm Character value specifying hic data normalization method. This value must be found in the .hic file. Default value is \code{norm = "KR"}.
#' @param matrix Character value indicating the type of matrix to output. Default value is \code{matrix = "observed"}. Options are:
#' \itemize{
#' \item{\code{"observed"}: }{Observed counts.}
#' \item{\code{"oe"}: }{Observed/expected counts.}
#' }
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#'
#' @return Returns a 3-column dataframe in sparse upper triangular format with the following columns: \code{chrom}, \code{altchrom}, \code{counts}.
#'
#' @seealso \link[strawr]{straw}
#'
#' @export
bb_readHic <- function(file, chrom, chromstart = NULL, chromend = NULL, altchrom = NULL, altchromstart = NULL, altchromend = NULL, assembly = "hg19", resolution = "auto", res_scale = "BP",
                       zrange = NULL, norm = "KR",  matrix = "observed", params = NULL){


  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_rhic
  errorcheck_bb_rhic <- function(hic, chrom, chromstart, chromend, zrange, altchrom, altchromstart, altchromend, norm, res_scale, assembly){

    ## hic input needs to be a path to a .hic file
    if (class(hic) != "character"){

      stop("Invalid input. Input needs to be a path to a .hic file.", call. = FALSE)

    }

    if ((file_ext(hic) != "hic")){

      stop("Invalid input. File must have a \".hic\" extension", call. = FALSE)

    }

    if (!file.exists(hic)){

      stop(paste("File", hic, "does not exist."), call. = FALSE)
    }


    ## Can't have only one NULL chromstart or chromend
    if ((is.null(chromstart) & !is.null(chromend)) | (is.null(chromend) & !is.null(chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.", call. = FALSE)

    }

    ## Not supporting chrM
    if (chrom == "chrM"){

      stop("chrM not supported.", call. = FALSE)

    }

    ## Even though straw technically works without "chr" for hg19, will not accept for consistency purposes
    if (assembly == "hg19"){

      if (grepl("chr", chrom) == FALSE){

        stop(paste(paste0("'",chrom, "'"), "is an invalid input for an hg19 chromsome. Please specify chromosome as", paste0("'chr", chrom, "'.")), call. = FALSE)
      }

    }


    if (!is.null(chromstart) & !is.null(chromend)){

      ## Chromstart should be smaller than chromend
      if (chromstart > chromend){

        stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)

      }

    }

    if (!is.null(altchrom)){

      if (altchrom == "chrM"){
        stop("chrM not supported.", call. = FALSE)

      }

      ## Even though straw technically works without "chr" for hg19, will not accept for consistency purposes
      if (assembly == "hg19"){

        if (grepl("chr", altchrom) == FALSE){

          stop(paste(paste0("'",altchrom, "'"), "is an invalid input for an hg19 chromsome. Please specify chromosome as", paste0("'chr", altchrom, "'.")), call. = FALSE)
        }

      }

      ## Can't specify altchrom without a chrom
      if (is.null(chrom)){

        stop("Specified \'altchrom\', but did not give \'chrom\'.", call. = FALSE)

      }

      ## Can't have only one NULL altchromstart or altchromend
      if ((is.null(altchromstart) & !is.null(altchromend)) | (is.null(altchromend) & !is.null(altchromstart))){

        stop("Cannot have one \'NULL\' \'altchromstart\' or \'altchromend\'.", call. = FALSE)

      }
      ## Altchromstart should be smaller than altchromend
      if (!is.null(altchromstart) & !is.null(altchromend)){
        if (altchromstart > altchromend){

          stop("\'altchromstart\' should not be larger than \'altchromend\'.", call. = FALSE)

        }

      }

      ## If giving same chrom and altchrom, need to specify chromstart/chromend and altchromstart/altchromend

      if (chrom == altchrom){

        if (is.null(chromstart) | is.null(chromend) | is.null(altchromstart) | is.null(altchromend)){

          stop("If giving the same \'chrom\' and \'altchrom\', please specify \'chromstart\', \'chromend\', \'altchromstart\', and \'altchromend\'.
               If trying to get all interactions between one chromosome, just specify \'chrom\'.", call. = FALSE)

        }

      }

    }

    ## Ensure properly formatted zrange
    if (!is.null(zrange)){

      ## zrange needs to be a vector
      if (!is.vector(zrange)){

        stop("\'zrange\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be length 2
      if (length(zrange) != 2){

        stop("\'zrange\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(zrange)){

        stop("\'zrange\' must be a vector of two numbers.", call. = FALSE)

      }

      ## second value should be larger than the first value
      if (zrange[1] >= zrange[2]){

        stop("\'zrange\' must be a vector of two numbers in which the 2nd value is larger than the 1st.", call. = FALSE)

      }

    }


    ## Check for valid "res_scale" parameter
    if(!(res_scale %in% c("BP", "FRAG"))){

      stop("Invalid \'res_scale\'.  Options are \'BP\' and \'FRAG\'.", call. = FALSE)

    }

  }

  ## Define a function that determines a best resolution for size of region
  auto_resolution <- function(chromstart, chromend){

    if (is.null(chromstart) & is.null(chromend)){
      autoRes <- 500000
    } else {
      dataRange <- chromend - chromstart
      if (dataRange >= 150000000){
        autoRes <- 500000
      } else if (dataRange >= 75000000 & dataRange < 150000000){
        autoRes <- 250000
      } else if (dataRange >= 35000000 & dataRange < 75000000){
        autoRes <- 100000
      } else if (dataRange >= 20000000 & dataRange < 35000000){
        autoRes <- 50000
      } else if (dataRange >= 5000000 & dataRange < 20000000){
        autoRes <- 25000
      } else if (dataRange >= 3000000 & dataRange < 5000000){
        autoRes <- 10000
      } else {
        autoRes <- 5000
      }
    }

    return(as.integer(autoRes))

  }

  ## Define a function to parse chromsome/region for Straw
  parse_region <- function(chrom, chromstart, chromend, assembly){

    if (assembly == "hg19"){
      strawChrom <- gsub("chr", "", chrom)
    } else {
      strawChrom <- chrom
    }

    if (is.null(chromstart) & is.null(chromend)){

      regionStraw <- strawChrom

    } else {

      ## Keep chromstart and chromend without scientific notation for processing with Straw

      chromstart <- format(chromstart, scientific = FALSE)
      chromend <- format(chromend, scientific = FALSE)
      regionStraw <- paste(strawChrom, chromstart, chromend, sep = ":")

    }

    return(regionStraw)

  }

  ## Define a function to reorder chromsomes to put "chrom" input in col1
  orderChroms <- function(hic, chrom, altchrom, assembly){

    if (assembly == "hg19"){
      chrom <- gsub("chr", "", chrom)
      altchrom <- gsub("chr", "", altchrom)
    }

    if (!"X" %in% chrom & !"Y" %in% chrom){
      chrom <- as.numeric(chrom)
    }

    if (!"X" %in% altchrom & !"Y" %in% altchrom){
      altchrom <- as.numeric(altchrom)
    }

    ## CASE 1: two numbers
    if (all(c(class(chrom), class(altchrom)) == "numeric")){

      if (chrom > altchrom){
        hic <- hic[, c(2, 1, 3)]
      }


    } else if (any(c(class(chrom), class(altchrom)) == "numeric")){
    ## CASE 2: number and X/Y
      if (class(altchrom) == "numeric"){
        hic <- hic[, c(2, 1, 3)]
      }

    } else {
    ## CASE 3: X and Y
      if ("Y" %in% chrom){
        hic <- hic[, c(2, 1, 3)]
      }

    }

    return(hic)

  }

  ## Define a function to scale data with zrange
  scale_data <- function(upper, zrange){

    if (!is.null(zrange)){

      upper$counts[upper$counts <= zrange[1]] <- zrange[1]
      upper$counts[upper$counts >= zrange[2]] <- zrange[2]

    } else {

      #if null, zrange will be set to (0, max(data))
      upper$counts[upper$counts <= 0] <- 0

    }

    return(upper)

  }

  ## Define a function to rename columns
  rename_columns <- function(upper, chrom, altchrom){

    if (is.null(altchrom)){

      colnames(upper) <- c(chrom, chrom, "counts")

    } else {

      colnames(upper) <- c(chrom, altchrom, "counts")

    }

    return(upper)
  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(resolution)) resolution <- NULL
  if(missing(norm)) norm <- NULL
  if(missing(res_scale)) res_scale <- NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(matrix)) matrix <- NULL

  ## Check if hic/chrom arguments are missing (could be in object)
  if(!hasArg(file)) file <- NULL
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_rhic <- structure(list(file = file, chrom = chrom, chromstart = chromstart, chromend = chromend, resolution = resolution, zrange = zrange,
                            norm = norm, res_scale = res_scale, assembly = assembly, matrix = matrix, altchrom = altchrom, altchromstart = altchromstart, altchromend = altchromend), class = "bb_rhic")

  bb_rhic <- parseParams(bb_params = params, object_params = bb_rhic)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_rhic$resolution)) bb_rhic$resolution <- "auto"
  if(is.null(bb_rhic$norm)) bb_rhic$norm <- "KR"
  if(is.null(bb_rhic$res_scale)) bb_rhic$res_scale <- "BP"
  if(is.null(bb_rhic$assembly)) bb_rhic$assembly <- "hg19"
  if(is.null(bb_rhic$matrix)) bb_rhic$matrix <- "observed"

  if(is.null(bb_rhic$file)) stop("argument \"file\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rhic$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  bb_rhic$assembly <- parse_bbAssembly(assembly = bb_rhic$assembly)

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_rhic(hic = bb_rhic$file, chrom = bb_rhic$chrom, chromstart = bb_rhic$chromstart, chromend = bb_rhic$chromend, zrange = bb_rhic$zrange, altchrom = bb_rhic$altchrom,
                     altchromstart = bb_rhic$altchromstart, altchromend = bb_rhic$altchromend, norm = bb_rhic$norm, res_scale = bb_rhic$res_scale, assembly = bb_rhic$assembly$Genome)

  # ======================================================================================================================================================================================
  # SET PARAMETERS
  # ======================================================================================================================================================================================

  parse_chromstart <- bb_rhic$chromstart
  parse_chromend <- bb_rhic$chromend
  parse_altchromstart <- bb_rhic$altchromstart
  parse_altchromend <- bb_rhic$altchromend


  ## For off diagonal plotting, grabbing whole symmetric region
  if (!is.null(bb_rhic$altchrom)){

    if (bb_rhic$chrom == bb_rhic$altchrom){

      parse_chromstart <- min(bb_rhic$chromstart, bb_rhic$altchromstart)
      parse_chromend <- max(bb_rhic$chromend, bb_rhic$altchromend)

    }

  }

  # ======================================================================================================================================================================================
  # PARSE REGIONS
  # ======================================================================================================================================================================================

  chromRegion <- parse_region(chrom = bb_rhic$chrom, chromstart = parse_chromstart, chromend = parse_chromend, assembly = bb_rhic$assembly$Genome)

  if (is.null(bb_rhic$altchrom)){

    altchromRegion <- chromRegion

  } else {

    if (bb_rhic$chrom == bb_rhic$altchrom){

      altchromRegion <- chromRegion

    } else {

      altchromRegion <- parse_region(chrom = bb_rhic$altchrom, chromstart = parse_altchromstart, chromend = parse_altchromend, assembly = bb_rhic$assembly$Genome)

    }

  }

  # ======================================================================================================================================================================================
  # ADJUST RESOLUTION
  # ======================================================================================================================================================================================

  if (bb_rhic$resolution == "auto"){
    bb_rhic$resolution <- auto_resolution(chromstart = bb_rhic$chromstart, chromend = bb_rhic$chromend)
    bb_rhic$res_scale <- "BP"
  }

  # ======================================================================================================================================================================================
  # EXTRACT SPARSE UPPER TRIANGULAR USING STRAW
  # ======================================================================================================================================================================================

  errorFunction <- function(c){
    upper <- data.frame(matrix(nrow = 0, ncol = 3))
    colnames(upper) <- c("x", "y", "counts")
    return(upper)
  }


  upper <-
    tryCatch(strawr::straw(bb_rhic$matrix, bb_rhic$norm, bb_rhic$file, toString(chromRegion), toString(altchromRegion), bb_rhic$res_scale, bb_rhic$resolution),
             error = errorFunction)

  # ======================================================================================================================================================================================
  # REORDER COLUMNS BASED ON CHROM/ALTCHROM INPUT
  # ======================================================================================================================================================================================

  if (!is.null(bb_rhic$altchrom)){

    upper <- orderChroms(hic = upper, chrom = bb_rhic$chrom, altchrom = bb_rhic$altchrom, assembly = bb_rhic$assembly$Genome)
    colnames(upper) <- c("x", "y", "counts")
  }

  # ======================================================================================================================================================================================
  # SCALE DATA WITH ZRANGE
  # ======================================================================================================================================================================================

  scaled_data <- scale_data(upper = upper, zrange = bb_rhic$zrange)

  # ======================================================================================================================================================================================
  # FORMAT DATA IN PROPER ORDER AND WITH LABELS
  # ======================================================================================================================================================================================

  renamed_data <- rename_columns(upper = scaled_data, chrom = bb_rhic$chrom, altchrom = bb_rhic$altchrom)

  # ======================================================================================================================================================================================
  # REMOVE NAN VALUES
  # ======================================================================================================================================================================================

  renamed_data <- na.omit(renamed_data)

  # ======================================================================================================================================================================================
  # RETURN DATAFRAME
  # ======================================================================================================================================================================================
  if (nrow(renamed_data) == 0){
    warning("No data found in region.", call. = FALSE)
  } else {
    message(paste("Read in hic file with", bb_rhic$norm, "normalization at", bb_rhic$resolution, bb_rhic$res_scale, "resolution."))
  }
  return(renamed_data)
}
