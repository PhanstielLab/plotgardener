#' plots HiC interaction matrix
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param half what sides of square plot; options are "both", top", or "bottom"
#' @param resolution the width in bp of each pixel; options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param just a string or numeric vector specifying the justification of the viewport relative to its (x, y) location
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#' @param draw A logical value indicating whether graphics output should be produced
#' @param altchrom alternate chromosome for off-diagonal plotting or interchromosomal plotting
#' @param altchromstart alternate chromosome start for off-diagonal plotting or interchromosomal plotting
#' @param altchromend alternate chromosome end for off-diagonal plotting or interchromosomal plotting
#' @param althalf if plotting altchrom region, which side off diagonal to plot; options are "top" or "bottom"
#' @param norm if giving .hic file, hic data normalization; must be found in hic file
#'
#' @return Function will plot a HiC interaction matrix and return a bb_hicPlot object
#'
#' @author Nicole Kramer
#'
#' @examples
#'
#' @export
#'
#'
bb_plotHic <- function(hic, chrom = 1, chromstart = 22500000, chromend = 23200000, half = "both", resolution = 10000, zrange = NULL,
                       palette = colorRampPalette(c("white", "dark red")), width = NULL, height = NULL, x = NULL, y = NULL,
                       just = c("left", "top"), default.units = "inches", draw = TRUE, altchrom = NULL, altchromstart = NULL, altchromend = NULL, althalf = NULL,
                       norm = "KR", ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_plothic
  errorcheck_bb_plothic <- function(hic, hic_plot, norm){

    ###### hic/norm #####

    ## if it's a dataframe or datatable, it needs to be properly formatted
    if ("data.frame" %in% class(hic) && ncol(hic) != 3){

      stop("Invalid dataframe format.  Input a dataframe with 3 columns: chrA, chrB, counts.", call. = FALSE)

    }

    if (!"data.frame" %in% class(hic)){

      ## if it's a file path, it needs to be a .hic file
      if (file_ext(hic) != "hic"){

        stop("Invalid input. File must have a \".hic\" extension", call. = FALSE)

      }

      ## if it's a file path, it needs to exist
      if (!file.exists(hic)){

        stop(paste("File", hic, "does not exist."), call. = FALSE)

      }

      ## if it's a valid .hic file, it needs to have a valid norm parameter
      if (is.null(norm)){

        stop("If providing .hic file, please specify \'norm\'.", call. = FALSE)

      }

    }

    ###### chrom/chromstart/chromend/altchrom/altchromstart/altchromend #####

    ## Can't have only one NULL chromstart or chromend
    if (any(is.null(hic_plot$chromstart), is.null(hic_plot$chromend))){

      stop("Please specify \'chromstart\' and \'chromend\'.", call. = FALSE)

    }

    ## Chromstart should be smaller than chromend
    if (hic_plot$chromstart > hic_plot$chromend){

      stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)

    }

    if (!is.null(hic_plot$altchrom)){

      ## Can't specify altchrom without a chrom
      if (is.null(hic_plot$chrom)){

        stop("Specified \'altchrom\', but did not give \'chrom\'.", call. = FALSE)

      }

      ## Can't have only one NULL altchromstart or altchromend
      if (any(is.null(hic_plot$altchromstart), is.null(hic_plot$altchromend))){

        stop("If specifying alternate chromosome, need to give \'altchromstart\' and \'altchromend\'.", call. = FALSE)

      }

      ## Altchromstart should be smaller than altchromend
      if (hic_plot$altchromstart > hic_plot$altchromend){

        stop("\'altchromstart\' should not be larger than \'altchromend\'.", call. = FALSE)

      }

      ## Check to see if region is square
      if ((hic_plot$chromend - hic_plot$chromstart) != (hic_plot$altchromend - hic_plot$altchromstart)){

        warning("Trying to plot non-square region.", call. = FALSE)

      }

    }

    ###### zrange #####

    ## Ensure properly formatted zrange
    if (!is.null(hic_plot$zrange)){

      ## zrange needs to be a vector
      if (!is.vector(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be length 2
      if (length(hic_plot$zrange) != 2){

        stop("\'zrange\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of two numbers.", call. = FALSE)

      }

      ## second value should be larger than the first value
      if (hic_plot$zrange[1] >= hic_plot$zrange[2]){

        stop("\'zrange\' must be a vector of two numbers in which the 2nd value is larger than the 1st.", call. = FALSE)

      }

    }

    ###### resolution #####

    if (is.null(hic_plot$resolution)){

      stop("Invalid \'resolution\' value.  Options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000.", call. = FALSE)

    } else {

      if (!(hic_plot$resolution %in% c(2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000))){

        stop("Invalid \'resolution\' value.  Options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000.", call. = FALSE)

      }

    }

    ###### half/althalf #####

    if (is.null(hic_plot$altchrom)){

      if (!(hic_plot$half %in% c("both", "top", "bottom"))){

        stop("Invalid \'half\'.  Options are \'both\', \top\', or \'bottom\'.", call. = FALSE)

      }

    } else {

      if (is.null(hic_plot$althalf)){

        stop("If specifying alternate chromosome, please provide \'althalf\'.", call. = FALSE)
      }

      if (!hic_plot$althalf %in% c("top", "bottom")){

        stop("Invalid \'althalf\'.  Options are \'top\' or \'bottom\'.", call. = FALSE)

      }

      if (!is.null(hic_plot$half)){

        warning("\'half\' option is unused.  \'althalf\' option will be used.", call. = FALSE)

      }


    }

  }

  ## Define a function to check range of data in dataframe
  check_dataframe <- function(hic, hic_plot){

    if (is.null(hic_plot$altchrom)){

      if (min(hic[,1]) > hic_plot$chromstart | max(hic[,1]) < hic_plot$chromend | min(hic[,2]) > hic_plot$chromstart | max(hic[,2]) < hic_plot$chromend){

        warning("Data is incomplete for the specified range.", call. = FALSE)

      }

    } else {

      if (min(hic[,1]) > hic_plot$chromstart | max(hic[,1]) < hic_plot$chromend | min(hic[,2]) > hic_plot$altchromstart | max(hic[,2]) < hic_plot$altchromend){

        warning("Data is incomplete for the specified range.", call. = FALSE)
      }

    }

  }

  ## Define a function that reads in hic data for bb_plothic
  read_data <- function(hic, hic_plot, norm){

    ## if .hic file, read in with bb_rhic
    if (!(class(hic) %in% "data.frame")){

      message(paste("Reading in hic file with", norm, "normalization."))

      readchromstart <- hic_plot$chromstart - hic_plot$resolution
      readchromend <- hic_plot$chromend + hic_plot$resolution
      readaltchromstart <- hic_plot$altchromstart - hic_plot$resolution
      readaltchromend <- hic_plot$altchromend + hic_plot$resolution

      hic <- bb_readHic(hic = hic, chrom = hic_plot$chrom, chromstart = readchromstart, chromend = readchromend,
                     resolution = hic_plot$resolution, zrange = hic_plot$zrange, norm = norm,
                     altchrom = hic_plot$altchrom, altchromstart = readaltchromstart,
                     altchromend = readaltchromend)

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

      hic <- hic[which(hic[,1] >= hic_plot$chromstart - hic_plot$resolution &
                         hic[,1] <= hic_plot$chromend + hic_plot$resolution &
                         hic[,2] >= hic_plot$chromstart - hic_plot$resolution &
                         hic[,2] <= hic_plot$chromend + hic_plot$resolution),]

    } else {

      if (hic_plot$chrom != hic_plot$altchrom){

        hic <- hic[which(hic[,1] >= hic_plot$chromstart - hic_plot$resolution &
                           hic[,1] <= hic_plot$chromend + hic_plot$resolution &
                           hic[,2] >= hic_plot$altchromstart - hic_plot$resolution &
                           hic[,2] <= hic_plot$altchromend + hic_plot$resolution),]

      }


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

      xscale <- c(hic_plot$chromstart, hic_plot$chromend)
      yscale <- xscale

    } else {

      if (hic_plot$althalf == "bottom"){

        xscale <- c(hic_plot$chromstart, hic_plot$chromend)
        yscale <- c(hic_plot$altchromstart, hic_plot$altchromend)

      }

      if (hic_plot$althalf == "top"){

        xscale <- c(hic_plot$altchromstart, hic_plot$altchromend)
        yscale <- c(hic_plot$chromstart, hic_plot$chromend)

      }


    }

    return(list(xscale, yscale))
  }

  ## Define a function that subsets the hic dataframe into which will be squares and which will be triangles
  hic_shapes <- function(hic, hic_plot){

    if (!is.null(hic_plot$altchrom)){

      ## all squares
      squares <- hic
      triangles <- NULL

    } else {

      if (hic_plot$half == "both"){

        ## all squares
        squares <- hic
        triangles <- NULL

      } else if (hic_plot$half == "top"){

        ## squares for top half
        ## triangles for diagonal

        squares <- hic[which(hic[,2] > hic[,1]),]
        triangles <- hic[which(hic[,2] == hic[,1]),]

      } else if (hic_plot$half == "bottom"){

        ## squares for bottom half
        ## triangles for diagonal

        squares <- hic[which(hic[,2] < hic[,1]),]
        triangles <- hic[which(hic[,2] == hic[,1]),]


      }

    }

  return(list(squares, triangles))

  }

  ## Define a function that makes grobs for the hic diagonal
  hic_diagonal <- function(hic, hic_plot){

    col <- hic[4]
    x <- as.numeric(hic[1])
    y <- as.numeric(hic[2])

    xleft = x
    xright = x + hic_plot$resolution
    ybottom = y
    ytop = y + hic_plot$resolution

    if (hic_plot$half == "top"){

      hic_triangle <- polygonGrob(x = c(xleft, xleft, xright),
                                  y = c(ybottom, ytop, ytop),
                                  gp = gpar(col = NA, fill = col),
                                  default.units = "native")


    } else if (hic_plot$half == "bottom"){


      hic_triangle <- polygonGrob(x = c(xleft, xright, xright),
                                  y = c(ybottom, ybottom, ytop),
                                  gp = gpar(col = NA, fill = col),
                                  default.units = "native")


      }

    assign("hic_grobs", addGrob(gTree = get("hic_grobs", envir = bbEnv), child = hic_triangle), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  hic_plot <- structure(list(chrom = chrom, chromstart = as.numeric(chromstart), chromend = as.numeric(chromend), altchrom = altchrom, altchromstart = altchromstart,
                             altchromend = altchromend, x = x, y = y, width = width, height = height, justification = just,
                             zrange = zrange, resolution = resolution, half = half, althalf = althalf, color_palette = NULL, grobs = NULL), class = "bb_hic")
  attr(x = hic_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = hic_plot)
  errorcheck_bb_plothic(hic = hic, hic_plot = hic_plot, norm = norm)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  hic_plot <- defaultUnits(object = hic_plot, default.units = default.units)

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  hic <- read_data(hic = hic, hic_plot = hic_plot, norm = norm)

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
    warning("Can't scale data to colors.", call. = FALSE)

  }


  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport xscale and yscale
  scale <- vp_scale(hic_plot = hic_plot)

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_hic", length(grep(pattern = "bb_hic", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(1, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = scale[[1]], yscale = scale[[2]],
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      vp$name <- "bb_hic1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = hic_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = scale[[1]], yscale = scale[[2]],
                   just = just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("hic_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  ## Determine which grobs will be squares [[1]] and which will be triangles [[2]]
  shapes <- hic_shapes(hic = hic, hic_plot = hic_plot)

  if (!is.null(shapes[[1]])){

    ## Make square grobs and add to grob gTree
    hic_squares <- rectGrob(x = shapes[[1]]$x,
                            y = shapes[[1]]$y,
                            just = c("left", "bottom"),
                            width = resolution,
                            height = resolution,
                            gp = gpar(col = NA, fill = shapes[[1]]$color),
                            default.units = "native")

    assign("hic_grobs", addGrob(gTree = get("hic_grobs", envir = bbEnv), child = hic_squares), envir = bbEnv)

  }

  ## Make triangle grobs and add to grob gTree
  if (!is.null(shapes[[2]])){

    invisible(apply(shapes[[2]], 1, hic_diagonal, hic_plot = hic_plot))

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("hic_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  hic_plot$grobs <- get("hic_grobs", envir = bbEnv)

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
