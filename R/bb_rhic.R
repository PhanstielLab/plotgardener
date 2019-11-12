#' extracts HiC data from .hic file using Straw
#'
#'
#' @param hic path to .hic file
#' @param chrom if not alternative chromosome, chromosome of desired region
#' @param chromstart chromosome start position of chrom
#' @param chromend chromosome end position of chrom
#' @param resolution the width of each pixel
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min; if null, zrange will be set to (0, max(data))
#' @param norm hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#' @param res_scale resolution scale; options are "BP" and "FRAG"
#' @param altchrom if looking at region between two different chromosomes, this is the specified alternative chromosome
#' @param altchromstart if looking at region between two different chromosomes, start position of altchrom
#' @param altchromend if looking at region between two different chromsomes, end position of altchrom
#'
#'
#' @export

bb_rhic <- function(hic, chrom, chromstart = NULL, chromend = NULL, resolution = 10000, zrange = NULL,
                    norm = "KR", res_scale = "BP", altchrom = NULL, altchromstart = NULL, altchromend = NULL){


  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_rhic
  errorcheck_bb_rhic <- function(hic, chromstart, chromend, zrange, altchrom, altchromstart, altchromend, norm, res_scale){

    ## hic input needs to be a path to a .hic file
    if (class(hic) != "character"){

      stop("Invalid input. Input needs to be a path to a .hic file.")

    }

    if ((file_ext(hic) != "hic")){

      stop("Invalid input. File must have a \".hic\" extension")

    }

    if (!file.exists(hic)){

      stop(paste("File", hic, "does not exist."))
    }

    ## Can't have only one NULL chromstart or chromend
    if ((is.null(chromstart) & !is.null(chromend)) | (is.null(chromend) & !is.null(chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.")

    }

    if (!is.null(chromstart) & !is.null(chromend)){

      ## Chromstart should be smaller than chromend
      if (chromstart > chromend){

        stop("\'chromstart\' should not be larger than \'chromend\'.")

      }

    }

    if (!is.null(altchrom)){

      ## Can't specify altchrom without a chrom
      if (is.null(chrom)){

        stop("Specified \'altchrom\', but did not give \'chrom\'.")

      }

      ## Can't have only one NULL altchromstart or altchromend
      if ((is.null(altchromstart) & !is.null(altchromend)) | (is.null(altchromend) & !is.null(altchromstart))){

        stop("Cannot have one \'NULL\' \'altchromstart\' or \'altchromend\'.")

      }
      ## Altchromstart should be smaller than altchromend

      if (altchromstart > altchromend){

        stop("\'altchromstart\' should not be larger than \'altchromend\'.")

      }

      ## If giving same chrom and altchrom, need to specify chromstart/chromend and altchromstart/altchromend

      if (chrom == altchrom){

        if (is.null(chromstart) | is.null(chromend) | is.null(altchromstart) | is.null(altchromend)){

          stop("If giving the same \'chrom\' and \'altchrom\', please specify \'chromstart\', \'chromend\', \'altchromstart\', and \'altchromend\'.
               If trying to get all interactions between one chromosome, just specify \'chrom\'.")

        }

      }

    }

    ## Ensure properly formatted zrange
    if (!is.null(zrange)){

      ## zrange needs to be a vector
      if (!is.vector(zrange)){

        stop("\'zrange\' must be a vector of length 2.")

      }

      ## zrange vector needs to be length 2
      if (length(zrange) != 2){

        stop("\'zrange\' must be a vector of length 2.")

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(zrange)){

        stop("\'zrange\' must be a vector of two numbers.")

      }

      ## second value should be larger than the first value
      if (zrange[1] >= zrange[2]){

        stop("\'zrange\' must be a vector of two numbers in which the 2nd value is larger than the 1st.")

      }

    }

    ## Check for valid "norm" parameter
    if (!(norm %in% c("NONE", "VC", "VC_SQRT", "KR"))){

      stop("Invalid \'norm\'.  Options are \'NONE\', \'VC\', \'VC_SQRT\', or \'KR\'.")

    }

    ## Check for valid "res_scale" parameter
    if(!(res_scale %in% c("BP", "FRAG"))){

      stop("Invalid \'res_scale\'.  Options are \'BP\' and \'FRAG\'.")

    }

  }

  ## Define a function to parse chromsome/region for Straw
  parse_region <- function(chrom, chromstart, chromend){

    if (is.null(chromstart) & is.null(chromend)){

      regionStraw <- chrom

    } else {

      ## Keep chromstart and chromend without scientific notation for processing with Straw

      chromstart <- format(chromstart, scientific = FALSE)
      chromend <- format(chromend, scientific = FALSE)
      regionStraw <- paste(chrom, chromstart, chromend, sep = ":")

    }

    return(regionStraw)

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
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_rhic(hic = hic, chromstart = chromstart, chromend = chromend, zrange = zrange, altchrom = altchrom, altchromstart = altchromstart,
                     altchromend = altchromend, norm = norm, res_scale = res_scale)

  # ======================================================================================================================================================================================
  # SET PARAMETERS
  # ======================================================================================================================================================================================
  parse_chromstart <- chromstart
  parse_chromend <- chromend
  parse_altchromstart <- altchromstart
  parse_altchromend <- altchromend


  ## For off diagonal plotting, grabbing whole symmetric region
  if (!is.null(altchrom)){

    if (chrom == altchrom){

      parse_chromstart <- min(chromstart, altchromstart)
      parse_chromend <- max(chromend, altchromend)

    }

  }

  # ======================================================================================================================================================================================
  # PARSE REGIONS
  # ======================================================================================================================================================================================

  chromRegion <- parse_region(chrom = chrom, chromstart = parse_chromstart, chromend = parse_chromend)

  if (is.null(altchrom)){

    altchromRegion <- chromRegion

  } else {

    if (chrom == altchrom){

      altchromRegion <- chromRegion

    } else {

      altchromRegion <- parse_region(chrom = altchrom, chromstart = parse_altchromstart, chromend = parse_altchromend)

    }

  }

  # ======================================================================================================================================================================================
  # EXTRACT SPARSE UPPER TRIANGULAR USING STRAW
  # ======================================================================================================================================================================================

  upper <- straw_R(sprintf("%s %s %s %s %s %i", norm, hic, chromRegion, altchromRegion, res_scale, resolution))

  # ======================================================================================================================================================================================
  # REORDER COLUMNS BASED ON CHROM/ALTCHROM INPUT
  # ======================================================================================================================================================================================
  if (!is.null(altchrom)){

    if (chrom > altchrom){

      upper <- upper[, c(2, 1, 3)]
      colnames(upper) <- c("x", "y", "counts")

    }

  }

  # ======================================================================================================================================================================================
  # SCALE DATA WITH ZRANGE
  # ======================================================================================================================================================================================

  scaled_data <- scale_data(upper = upper, zrange = zrange)

  # ======================================================================================================================================================================================
  # FILL IN MISSING GAPS OF DATA
  # ======================================================================================================================================================================================

  # complete_data <- fill_missing_data(dataframe = scaled_data, chrom = chrom, chromstart = chromstart, chromend = chromend, altchrom = altchrom,
  #                              altchromstart = altchromstart, altchromend = altchromend, resolution = resolution, fill_missing = fill_missing)

  # ======================================================================================================================================================================================
  # FORMAT DATA IN PROPER ORDER AND WITH LABELS
  # ======================================================================================================================================================================================

  #renamed_data <- rename_columns(upper = complete_data, chrom = chrom, altchrom = altchrom)
  renamed_data <- rename_columns(upper = scaled_data, chrom = chrom, altchrom = altchrom)

  # ======================================================================================================================================================================================
  # RETURN DATAFRAME
  # ======================================================================================================================================================================================
  if (nrow(renamed_data) == 0){

    warning("Warning: no data found in region.  Suggestions: check chromosome, check region.")
  }


  return(renamed_data)
}
