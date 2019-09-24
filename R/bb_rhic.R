#' extracts HiC data from .hic file using Straw
#'
#'
#' @param hic path to .hic file
#' @param chrom if not alternative chromosome, chromosome of desired region
#' @param chromstart chromosome start position of chrom, in bp
#' @param chromend chromosome end position of chrom, in bp
#' @param resolution the width in bp of each pixel
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param norm hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#' @param res_scale scale of normalization; options are "BP" and "FRAG"
#' @param altchrom if looking at region between two different chromosomes, this is the specified alternative chromsome
#' @param altchromstart if looking at region between two different chromosomes, start position of altchrom
#' @param altchromend if looking at region between two different chromsomes, end position of altchrom
#'
#'
#' @export

bb_rhic <- function(hic, chrom, chromstart = NULL, chromend = NULL, resolution = 10000, zrange = NULL,
                    norm = "KR", res_scale = "BP", altchrom = NULL, altchromstart = NULL, altchromend = NULL){

  # Parse chromosome and region in format for Straw
  # ======================================================================================================================================================================================
  if ((is.null(chromstart) & !is.null(chromend)) | (is.null(chromend) & !is.null(chromstart))){
    stop("Cannot have one \'NULL\' chromstart or chromend.")
  } else if (is.null(chromstart) & is.null(chromend)){
    regionChrom <- gsub(pattern = "chr", replacement = "", x = chrom)
    regionStraw <- regionChrom
  } else {
    regionChrom <- gsub(pattern = "chr", replacement = "", x = chrom)

    ## Keep chromstart and chromend without scientific notation for processing with Straw
    chromstart <- format(chromstart, scientific = FALSE)
    chromend <- format(chromend, scientific = FALSE)
    regionStraw <- paste(regionChrom, chromstart, chromend, sep = ":")
  }

  # Extract upper triangular using straw, depending on one chromsome interaction or multiple chromosome interactions
  # ======================================================================================================================================================================================
  if(is.null(altchrom)){

    upper <- straw_R(sprintf("%s %s %s %s %s %i", norm, hic, regionStraw, regionStraw, res_scale, resolution))

    ## Fill in gaps and complete data
    upper <- fill_missing_data(dataframe = upper, chrom = chrom, chromstart = chromstart, chromend = chromend, resolution = resolution)

    ## Rename column names with proper chromosome
    colnames(upper)[colnames(upper) == "x"] <- paste0("chr", regionChrom)
    colnames(upper)[colnames(upper) == "y"] <- paste0("chr", regionChrom)

  } else {
    if ((is.null(altchromstart) & !is.null(altchromend)) | (is.null(altchromend) & !is.null(altchromstart))){
      stop("Cannot have one \'NULL\' altchromstart or altchromend.")
    } else if (is.null(altchromstart) & is.null(altchromend)){
      regionChrom2 <- gsub(pattern = "chr", replacement = "", x = altchrom)
      regionStraw2 <- regionChrom2
    } else {
      regionChrom2 <- gsub(pattern = "chr", replacement = "", x = altchrom)

      ## Keep altchromstart and altchromend without scientific notation for processing with Straw
      altchromstart <- format(altchromstart, scientific = FALSE)
      altchromend <- format(altchromend, scientific = FALSE)
      regionStraw2 <- paste(regionChrom2, altchromstart, altchromend, sep = ":" )
    }

    upper <- straw_R(sprintf("%s %s %s %s %s %i", norm, hic, regionStraw, regionStraw2, res_scale, resolution))

    ## Fill in gaps and complete data
    upper <- fill_missing_data(dataframe = upper, chrom = chrom, chromstart = chromstart, chromend = chromend, altchrom = altchrom,
                                 altchromstart = altchromstart, altchromend = altchromend, resolution = resolution)

    ## Rename column names with proper chromosomes
    colnames(upper)[colnames(upper) == "x"] <- paste0("chr", min(as.numeric(regionChrom), as.numeric(regionChrom2)))
    colnames(upper)[colnames(upper) == "y"] <- paste0("chr", max(as.numeric(regionChrom), as.numeric(regionChrom2)))

  }


  ## Scale lower and upper bounds using zrange
  if(is.null(zrange)){

    zrange <- c(min(upper$counts), max(upper$counts))
    upper$counts[upper$counts <= zrange[1]] <- zrange[1]
    upper$counts[upper$counts >= zrange[2]] <- zrange[2]

  } else {
    stopifnot(is.vector(zrange), length(zrange) == 2, zrange[2] > zrange[1])
    upper$counts[upper$counts <= zrange[1]] <- zrange[1]
    upper$counts[upper$counts >= zrange[2]] <- zrange[2]
  }

  if (nrow(upper) == 0){
    warning("Warning: no data found in region.  Suggestions: check chromosome, check region")
  }

  return(as.data.frame(upper))
}
