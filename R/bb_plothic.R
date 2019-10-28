#' plots HiC interaction matrix
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param half what sides of square plot; options are "both", top", or "bottom"
#' @param resolution the width in bp of each pixel
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param width width of plot in inches
#' @param height height of plot in inches
#' @param x x-coordinate of where to place plot relative to top left of plot
#' @param y y-coordinate of where to place plot relative to top left of plot
#' @param just a string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param units units of height, width, and x and y location of the plot
#' @param altchrom alternate chromosome for off-diagonal plotting or interchromosomal plotting
#' @param altchromstart alternate chromosome start for off-diagonal plotting or interchromosomal plotting
#' @param altchromend alternate chromosome end for off-diagonal plotting or interchromosomal plotting
#' @param althalf if plotting altchrom region, which side off diagonal to plot; options are "top" or "bottom"
#' @param norm if giving .hic file, hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#' @param fill_missing option to fill missing gaps of data with 0's (TRUE) or NA's (FALSE)
#'
#' @details macOS Preview anti-aliases images and will make rasterized plot appear blurry.
#'
#' @return Function will plot a HiC interaction matrix and return a list of the zrange min, zrange max, and sequential color palette
#'
#' @author Nicole Kramer
#'
#' @examples
#'
#' @export
#'
#'
bb_plothic <- function(hic, chrom = 8, chromstart = 133600000, chromend = 134800000, half = "both", resolution = 10000, zrange = NULL,
                       palette = colorRampPalette(c("white", "dark red")), width = 3, height = 3, x = 1, y = 1,
                       just = c("left", "top"), units = "inches", altchrom = NULL, altchromstart = NULL, altchromend = NULL, althalf = NULL,
                       norm = "KR", fill_missing = T, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_plothic
  errorcheck_bb_plothic <- function(hic, hic_plot){

    # hic input
    # ===================================================================================
    ## If hic is data frame ensure that it is formatted correctly
    if (class(hic) %in% "data.frame" && ncol(hic) != 3){

      stop("Invalid dataframe format.  Input a dataframe with 3 columns: chrA, chrB, counts.")

    }

    ## If hic is path, make sure it has the correct extension and that it exists
    if (!class(hic) %in% "data.frame"){

      if (file_ext(hic) != "hic"){

        stop("Invalid input. File must have a \".hic\" extension")

      }

      if (!file.exists(hic)){

        stop(paste("File", hic, "does not exist."))

      }

    }

    # chrom/altchrom/chromstart/chromend/altchromstart/altchromend
    # ===================================================================================
    ## Can't have only one NULL chromstart or chromend
    if ((is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)) | (is.null(hic_plot$chromend) & !is.null(hic_plot$chromstart))){

      stop("Please specify \'chromstart\' and \'chromend\'.")

    }

    ## Chromstart should be smaller than chromend
    if (hic_plot$chromstart > hic_plot$chromend){

      stop("\'chromstart\' should not be larger than \'chromend\'.")

    }

    if (!is.null(hic_plot$altchrom)){

      ## Can't specify altchrom without a chrom
      if (is.null(hic_plot$chrom)){

        stop("Specified \'altchrom\', but did not give \'chrom\'.")

      }

      ## Can't have only one NULL altchromstart or altchromend
      if ((is.null(hic_plot$altchromstart) & !is.null(hic_plot$altchromend)) | (is.null(hic_plot$altchromend) & !is.null(hic_plot$altchromstart))){

        stop("If specifying alternate chromosome, need to give \'altchromstart\' and \'altchromend\'.")

      }

      ## Altchromstart should be smaller than altchromend
      if (hic_plot$altchromstart > hic_plot$altchromend){

        stop("\'altchromstart\' should not be larger than \'altchromend\'.")

      }

      ## Check to see if region is square
      if ((hic_plot$chromend - hic_plot$chromstart) != (hic_plot$altchromend - hic_plot$altchromstart)){

        warning("Trying to plot non-square region.")

      }

    }

    # zrange
    # ===================================================================================
    ## Ensure properly formatted zrange
    if (!is.null(hic_plot$zrange)){

      ## zrange needs to be a vector
      if (!is.vector(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of length 2.")

      }

      ## zrange vector needs to be length 2
      if (length(hic_plot$zrange) != 2){

        stop("\'zrange\' must be a vector of length 2.")

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of two numbers.")

      }

      ## second value should be larger than the first value
      if (hic_plot$zrange[1] >= hic_plot$zrange[2]){

        stop("\'zrange\' must be a vector of two numbers in which the 2nd value is larger than the 1st.")

      }

    }

    # parameter values
    # ===================================================================================
    if (class(hic_plot$additional_parameters$palette) != "function"){

      stop("Invalid palette. Please give a palette function.")

    }


    if (!(hic_plot$additional_parameters$half %in% c("both", "top", "bottom"))){

      stop("Invalid \'half\'.  Options are \'both\', \top\', or \'bottom\'.")

    }

    if (!is.null(hic_plot$altchrom) & is.null(hic_plot$additional_parameters$althalf)){

      stop("If specifying alternate chromosome, please provide \'althalf\'.")

    }

    if (!is.null(hic_plot$altchrom)){

      if (!(hic_plot$additional_parameters$althalf %in% c("top", "bottom"))){

        stop("Please specify \'althalf\'.  Options are \'top\' or \'bottom\'.")

      }

    }


    if (!class(hic) %in% "data.frame" & !(hic_plot$additional_parameters$norm %in% c("NONE", "VC", "VC_SQRT", "KR"))){

      stop("If providing .hic file, please specify \'norm\'.  Options are \'NONE\', \'VC\', \'VC_SQRT\', or \'KR\'.")

    }

  }

  ## Define a function to check range of data in dataframe
  check_dataframe <- function(hic, hic_plot){

    if (is.null(hic_plot$altchrom)){

      if (min(hic[,1]) > hic_plot$chromstart | max(hic[,1]) < hic_plot$chromend | min(hic[,2]) > hic_plot$chromstart | max(hic[,2]) < hic_plot$chromend){

        warning("Data is incomplete for the specified range.")

      }

    } else {

      if (min(hic[,1]) > hic_plot$chromstart | max(hic[,1]) < hic_plot$chromend | min(hic[,2]) > hic_plot$altchromstart | max(hic[,2]) < hic_plot$altchromend){

        warning("Data is incomplete for the specified range.")
      }

    }

  }

  ## Define a function that reads in hic data for bb_plothic
  read_data <- function(hic, hic_plot){

    ## if .hic file, read in with bb_rhic
    if (!(class(hic) %in% "data.frame")){

      message(paste("Reading in hic file with", norm, "normalization and fill_missing =", fill_missing))

      readchromstart <- hic_plot$chromstart - hic_plot$additional_parameters$resolution
      readchromend <- hic_plot$chromend + hic_plot$additional_parameters$resolution
      readaltchromstart <- hic_plot$altchromstart - hic_plot$additional_parameters$resolution
      readaltchromend <- hic_plot$altchromend + hic_plot$additional_parameters$resolution

      hic <- bb_rhic(hic = hic, chrom = hic_plot$chrom, chromstart = readchromstart, chromend = readchromend,
                     resolution = hic_plot$additional_parameters$resolution, zrange = hic_plot$zrange,norm = hic_plot$additional_parameters$norm,
                     altchrom = hic_plot$altchrom, altchromstart = readaltchromstart,
                     altchromend = readaltchromend, fill_missing = hic_plot$additional_parameters$fill_missing)

    } else {

      message("Reading in dataframe.  Assuming \'chrom\' in column1 and \'altchrom\' in column2.")

      ## check range of data in dataframe
      check_dataframe(hic = hic, hic_plot = hic_plot)

    }

    ## Rename columns for later processing
    colnames(hic) <- c("x", "y", "counts")

    return(hic)

  }

  ## Define a function that subsets data
  subset_data <- function(hic, hic_plot){

    if(is.null(hic_plot$altchrom)){

      hic <- hic[which(hic[,1] >= hic_plot$chromstart - hic_plot$additional_parameters$resolution &
                         hic[,1] <= hic_plot$chromend + hic_plot$additional_parameters$resolution &
                         hic[,2] >= hic_plot$chromstart - hic_plot$additional_parameters$resolution &
                         hic[,2] <= hic_plot$chromend + hic_plot$additional_parameters$resolution),]

    } else {

      hic <- hic[which(hic[,1] >= hic_plot$chromstart - hic_plot$additional_parameters$resolution &
                         hic[,1] <= hic_plot$chromend + hic_plot$additional_parameters$resolution &
                         hic[,2] >= hic_plot$altchromstart - hic_plot$additional_parameters$resolution &
                         hic[,2] <= hic_plot$altchromend + hic_plot$additional_parameters$resolution),]
    }


    return(hic)
  }

  ## Define a function that sets the zrange
  set_zrange <- function(hic, hic_plot){

    ## no zrange, only one value
    if (is.null(hic_plot$zrange) & length(unique(hic$counts)) == 1){

      zrange <- c(unique(hic$counts), unique(hic$counts))
      hic_plot$zrange <- zrange

    }

    ## no zrange, multiple values
    if (is.null(hic_plot$zrange) & length(unique(hic$counts)) > 1){

      zrange <- c(0, max(hic$counts))
      hic_plot$zrange <- zrange

    }

    return(hic_plot)

  }

  ## Define a function that sets viewport xscale and yscale
  vp_scale <- function(hic_plot){

    if(is.null(hic_plot$altchrom)){

      xscale <- c(hic_plot$chromstart, (hic_plot$chromend + hic_plot$additional_parameters$resolution))
      yscale <- xscale

    } else {

      if (hic_plot$additional_parameters$althalf == "bottom"){

        xscale <- c(hic_plot$chromstart, (hic_plot$chromend + hic_plot$additional_parameters$resolution))
        yscale <- c(hic_plot$altchromstart, (hic_plot$altchromend + hic_plot$additional_parameters$resolution))

      }

      if (hic_plot$additional_parameters$althalf == "top"){

        xscale <- c(hic_plot$altchromstart, (hic_plot$altchromend + hic_plot$additional_parameters$resolution))
        yscale <- c(hic_plot$chromstart, (hic_plot$chromend + hic_plot$additional_parameters$resolution))

      }


    }

    return(list(xscale, yscale))
  }

  ## Define a function that plots squares and triangles on hic plot
  drawpoly <- function(df, hic_plot){

    col = df[4]
    x = as.numeric(df[1])
    y = as.numeric(df[2])

    xleft = x
    xright = x + hic_plot$additional_parameters$resolution
    ybottom = y
    ytop = y + hic_plot$additional_parameters$resolution

    if (!is.null(hic_plot$altchrom)){


      grid.polygon(x = c(xleft, xleft, xright, xright),
                   y = c(ybottom, ytop, ytop, ybottom), gp = gpar(col = NA, fill = col), default.units = "native")

    } else {

      if (hic_plot$additional_parameters$half == "both"){

        ## Plot all squares
        grid.polygon(x = c(xleft, xleft, xright, xright),
                     y = c(ybottom, ytop, ytop, ybottom), gp = gpar(col = NA, fill = col), default.units = "native")
      } else if (hic_plot$additional_parameters$half == "top"){

        ## Plot triangles along diagonal and squares above
        if (y > x){
          grid.polygon(x = c(xleft, xleft, xright, xright),
                       y = c(ybottom, ytop, ytop, ybottom), gp = gpar(col = NA, fill = col), default.units = "native")
        } else if (y == x) {
          grid.polygon(x = c(xleft, xleft, xright),
                       y = c(ybottom, ytop, ytop), gp = gpar(col = NA, fill = col), default.units = "native")
        }

      } else if (hic_plot$additional_parameters$half == "bottom"){
        ## Plot triangles along diagonal and squares below
        if (y < x){
          grid.polygon(x = c(xleft, xleft, xright, xright),
                       y = c(ybottom, ytop, ytop, ybottom), gp = gpar(col = NA, fill = col), default.units = "native")
        } else if (y == x) {
          grid.polygon(x = c(xleft, xright, xright),
                       y = c(ybottom, ybottom, ytop), gp = gpar(col = NA, fill = col), default.units = "native")
        }
      }

    }

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  hic_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, x = x, y = y, height = height,
                             width = width, units = units, altchrom = altchrom, altchromstart = altchromstart,
                             altchromend = altchromend, zrange = zrange, color_palette = NULL, just = just,
                             additional_parameters = list(half = half,
                                                          resolution = resolution, palette = palette, althalf = althalf,
                                                          norm = norm, fill_missing = fill_missing)), class = "hic_plot")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plothic(hic = hic, hic_plot = hic_plot)

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  hic <- read_data(hic = hic, hic_plot = hic_plot)

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  hic <- subset_data(hic = hic, hic_plot = hic_plot)

  # ======================================================================================================================================================================================
  # MAKE SYMMETRIC
  # ======================================================================================================================================================================================

  hicFlip = hic[, c(2, 1, 3)]
  colnames(hicFlip) <- c("x", "y", "counts")
  hic <- unique(rbind(hic, hicFlip))
  colnames(hic) = c("x", "y", "counts")


  # ======================================================================================================================================================================================
  # SET ZRANGE AND SCALE DATA
  # ======================================================================================================================================================================================

  hic_plot <- set_zrange(hic = hic, hic_plot = hic_plot)
  hic$counts[hic$counts <= hic_plot$zrange[1]] <- hic_plot$zrange[1]
  hic$counts[hic$counts >= hic_plot$zrange[2]] <- hic_plot$zrange[2]

  # ======================================================================================================================================================================================
  # CONVERT NUMBERS TO COLORS
  # ======================================================================================================================================================================================

  ## if we don't have an appropriate zrange (even after setting it based on a null zrange), can't scale to colors
  if (!is.null(hic_plot$zrange) & length(unique(hic_plot$zrange)) == 2){

    hic$color <- bb_maptocolors(hic$counts, col = palette, num = 100, range = hic_plot$zrange)
    sorted_colors <- unique(hic[order(hic$counts),]$color)
    hic_plot$color_palette <- sorted_colors

  } else {

    ## If we still have a null zrange or a length(unique(zrange)) == 1, means we couldn't do it in setzrange above (empty data or data with only 1 value)
    warning("Can't scale data to colors.")

  }


  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = hic_plot)

  ## Get viewport xscale and yscale
  scale <- vp_scale(hic_plot = hic_plot)

  ## Make viewport
  vp <- viewport(height = unit(page_coords[[1]]$height, page_coords[[3]]), width = unit(page_coords[[1]]$width, page_coords[[3]]),
                 x = unit(page_coords[[1]]$x, page_coords[[3]]), y = unit((page_coords[[2]]-page_coords[[1]]$y), page_coords[[3]]),
                 clip = "on", xscale = scale[[1]], yscale = scale[[2]], just = just)

  pushViewport(vp)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  invisible(apply(hic, 1, drawpoly, hic_plot = hic_plot))

  ## Go back up a viewport
  upViewport()

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================
  if (is.null(altchrom)){

    hic_plot$altchrom = chrom
    hic_plot$altchromstart = chromstart
    hic_plot$altchromend = chromend

  }

  return(hic_plot)

  }
