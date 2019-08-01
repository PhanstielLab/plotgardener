#' extracts HiC data from .hic file using Straw
#'
#'
#' @param hic path to .hic file
#' @param chrom if not alternative chromosome, chromosome of desired region
#' @param chromstart chromosome start position of chrom, in bp
#' @param chromend chromosome end position of chrom, in bp
#' @param resolution the width in bp of each pixel
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param format format of data wanted, which will depend on what kind of plot is ultimately desired; options are "sparse" and "full"
#' @param norm hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#' @param res_scale scale of normalization; options are "BP" and "FRAG"
#' @param altchrom if looking at region between two different chromosomes, this is the specified alternative chromsome
#' @param altchromstart if looking at region between two different chromosomes, start position of altchrom
#' @param altchromend if looking at region between two different chromsomes, end position of altchrom
#'
#'
#' @export

bb_rhic <- function(hic, chrom, chromstart = NULL, chromend = NULL, resolution = 10000, zrange = NULL, format = "sparse",
                    norm = "KR", res_scale = "BP", altchrom = NULL, altchromstart = NULL, altchromend = NULL){

  # Parse chromosome and region in format for Straw
  # ======================================================================================================================================================================================
  if ((is.null(chromstart) & !is.null(chromend)) | (is.null(chromend) & !is.null(chromstart))){
    stop("Cannot have one \'NULL\' chromstart or chromend.")
  } else if (is.null(chromstart) & is.null(chromend)){
    regionChrom <- gsub(pattern = "chr|chrom|CHR|CHROM", replacement = "", x = chrom)
    regionStraw <- regionChrom
  } else {
    regionChrom <- gsub(pattern = "chr|chrom|CHR|CHROM", replacement = "", x = chrom)

    ## Keep chromstart and chromend without scientific notation for processing with Straw
    chromstart <- format(chromstart, scientific = FALSE)
    chromend <- format(chromend, scientific = FALSE)
    regionStraw <- paste(regionChrom, chromstart, chromend, sep = ":")
  }

  # Extract upper triangular using straw, depending on one chromsome interaction or multiple chromosome interactions
  # ======================================================================================================================================================================================
  if(is.null(altchrom)){
    upper <- straw_R(sprintf("%s %s %s %s %s %i", norm, hic, regionStraw, regionStraw, res_scale, resolution))

    # Rename "x" and "y" with chromosome
    colnames(upper)[colnames(upper) == "x"] <- paste0("chr", regionChrom)
    colnames(upper)[colnames(upper) == "y"] <- paste0("chr" ,regionChrom)
  } else {
    if ((is.null(altchromstart) & !is.null(altchromend)) | (is.null(altchromend) & !is.null(altchromstart))){
      stop("Cannot have one \'NULL\' altchromstart or altchromend.")
    } else if (is.null(altchromstart) & is.null(altchromend)){
      regionChrom2 <- gsub(pattern = "chr|chrom|CHR|CHROM", replacement = "", x = altchrom)
      regionStraw2 <- regionChrom2
    } else {
      regionChrom2 <- gsub(pattern = "chr|chrom|CHR|CHROM", replacement = "", x = altchrom)

      ## Keep altchromstart and altchromend without scientific notation for processing with Straw
      altchromstart <- format(altchromstart, scientific = FALSE)
      altchromend <- format(altchromend, scientific = FALSE)
      regionStraw2 <- paste(regionChrom2, altchromstart, altchromend, sep = ":" )
    }

    upper <- straw_R(sprintf("%s %s %s %s %s %i", norm, hic, regionStraw, regionStraw2, res_scale, resolution))
    # Rename "x" and "y" with chromosomes
    colnames(upper)[colnames(upper) == "x"] <- paste0("chr", min(as.numeric(regionChrom), as.numeric(regionChrom2)))
    colnames(upper)[colnames(upper) == "y"] <- paste0("chr", max(as.numeric(regionChrom), as.numeric(regionChrom2)))
  }

  # Full format: get symmetric data, complete missing values, replace NA's with 0's
  # ======================================================================================================================================================================================
  if(format == "full"){
    ## "complete" function cannot have duplicate column names, so temporarily rename columns back to x and y
    colnames(upper) <- c("x", "y", "counts")
    lower <- upper[ ,c(2,1,3)]
    colnames(lower) <- c("x", "y", "counts")
    combined <- unique(rbind(upper, lower))
    combinedComplete <- tidyr::complete(combined, x, y)
    combinedComplete$counts[is.na(combinedComplete$counts)] <- 0

    ## Rename "x" and "y" columns with chromosome(s) again
    if(is.null(altchrom)){
      colnames(combinedComplete)[colnames(combinedComplete) == "x"] <- paste0("chr", regionChrom)
      colnames(combinedComplete)[colnames(combinedComplete) == "y"] <- paste0("chr" ,regionChrom)
    } else {
      colnames(combinedComplete)[colnames(combinedComplete) == "x"] <- paste0("chr", min(as.numeric(regionChrom), as.numeric(regionChrom2)))
      colnames(combinedComplete)[colnames(combinedComplete) == "y"] <- paste0("chr", max(as.numeric(regionChrom), as.numeric(regionChrom2)))
    }

    # Scale lower and upper bounds using zrange
    if(is.null(zrange)){
        zrange <- c(min(combinedComplete$counts), max(combinedComplete$counts))
        combinedComplete$counts[combinedComplete$counts <= zrange[1]] <- zrange[1]
        combinedComplete$counts[combinedComplete$counts >= zrange[2]] <- zrange[2]

    } else {
      stopifnot(is.vector(zrange), length(zrange) == 2, zrange[2] > zrange[1])
      combinedComplete$counts[combinedComplete$counts <= zrange[1]] <- zrange[1]
      combinedComplete$counts[combinedComplete$counts >= zrange[2]] <- zrange[2]
    }

    if (nrow(combinedComplete) == 0){
      warning("Warning: no data found in region.  Suggestions: check chromosome, check region")
    }

    return(as.data.frame(combinedComplete))
  }

  # Sparse format: upper sparse triangular format, scale with new zrange if given
  # ======================================================================================================================================================================================
  if (format == "sparse"){
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
}
