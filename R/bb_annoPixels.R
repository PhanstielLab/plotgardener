#' Annotate pixels in a Hi-C plot
#'
#' @param plot Hi-C plot object from \code{bb_plotHicSquare} or
#' \code{bb_plotHicTriangle} on which to annotate pixels.
#' @param data A string specifying the BEDPE file path or a dataframe in BEDPE
#' format specifying pixel positions.
#' @param type Character value specifying type of annotation.
#' Default value is \code{type = "box"}. Options are:
#' \itemize{
#' \item{\code{"box"}: }{Boxes are drawn around each pixel.}
#' \item{\code{"circle"}: }{Circles are drawn around each pixel.}
#' \item{\code{"arrow"}: }{Arrows are drawn pointing to each pixel.}
#' }
#' @param half Character value specifying which half of hic plots
#' to annotate. Triangle Hi-C plots will always default to the entirety of
#' the triangular plot. Default value is \code{half = "inherit"}. Options are:
#' \itemize{
#' \item{\code{"inherit"}: }{Pixels will be annotated on the \code{half}
#' inherited by the input Hi-C plot.}
#' \item{\code{"both"}: }{Pixels will be annotated on both halves of the
#' diagonal of a square Hi-C plot.}
#' \item{\code{"top"}: }{Pixels will be annotated on the upper diagonal
#' half of a square Hi-C plot.}
#' \item{\code{"bottom"}: }{Pixels will be annotated ont the bottom diagonal
#' half of a square Hi-C plot.}
#' }
#' @param shift Numeric specifying the number of pixels on either end of
#' main pixel in a box or circle. Numeric specifying number of pixels
#' for the length of an arrow.
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#' @param quiet A logical indicating whether or not to print messages.
#'
#' @return Returns a \code{bb_pixel} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data and BEDPE data
#' data("bb_imrHicData")
#' data("bb_bedpeData")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 4.5, height = 4, default.units = "inches")
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- bb_plotHicSquare(data = bb_imrHicData, resolution = 10000,
#'                             zrange = c(0, 70),
#'                             chrom = "chr21",
#'                             chromstart = 28000000, chromend = 30300000,
#'                             x = 0.5, y = 0.5, width = 3, height = 3,
#'                             just = c("left", "top"),
#'                             default.units = "inches")
#'
#' ## Annotate loops of both sides of Hi-C plot with squares
#' pixels <- bb_annoPixels(plot = hicPlot, data = bb_bedpeData, type = "box",
#'                         half = "both")
#'
#' ## Annotate loops on one side of Hi-C plot with arrows
#' ## and the other side with circles
#' bb_pagePlotRemove(plot = pixels)
#' pixels1 <- bb_annoPixels(plot = hicPlot, data = bb_bedpeData,
#'                          type = "arrow", half = "top", shift = 8)
#' pixels2 <- bb_annoPixels(plot = hicPlot, data = bb_bedpeData,
#'                          type = "circle", half = "bottom")
#'
#' ## Annotate heatmap legend
#' bb_annoHeatmapLegend(plot = hicPlot,
#'                      x = 3.6, y = 0.5, width = 0.12, height = 1.2,
#'                      just = c("left", "top"), default.units = "inches")
#'
#' ## Annotate genome label
#' bb_annoGenomeLabel(plot = hicPlot, x = 0.5, y = 3.53, scale = "Mb",
#'                    just = c("left", "top"))
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @export
bb_annoPixels <- function(plot, data, type = "box", half = "inherit",
                          shift = 4, params = NULL, quiet = FALSE,...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================
  ## For more accurate calculation of sqrt(2)
  two <- mpfr(2, 120)

  ## Define a function to catch errors for bb_annoPixels
  errorcheck_bb_annoLoops <- function(hic, loops, half, type, quiet){

    ###### hic #####

    ## check type of input for hic
    if (!class(hic) %in% c("bb_hicSquare", "bb_hicTriangle", "bb_hicRectangle")){

      stop("Input plot must be a plot of class \'bb_hicSquare\', \'bb_hicTriangle\', or \'bb_hicRectangle\'.", call. = FALSE)

    }

    ###### loops #####

    ## if loops is a dataframe or datatable, it needs to be properly formatted
    if ("data.frame" %in% class(loops) && ncol(loops) < 6){

      stop("Invalid dataframe format. Dataframe must be in BEDPE format.", call. = FALSE)

    }

    if ("data.frame" %in% class(loops) && nrow(loops) < 1){

      stop("\'data\' input contains no values.", call. = FALSE)

    }


    ## if it's a file path, it needs to exist
    if (!"data.frame" %in% class(loops)){

      ## File existence
      if (!file.exists(loops)){

        stop(paste("File", loops, "does not exist."), call. = FALSE)

      }

    }

    ###### half #####

    ## half needs to be a valid option
    if (!half %in% c("inherit", "both", "top", "bottom")){

      stop("Invalid \'half\'.  Options are \'inherit\', \'both\', \'top\', or \'bottom\'.", call. = FALSE)

    }

    ## half needs to be able to align with what kind of hic plot is plotted
    if (class(hic) == "bb_hicSquare"){

      if (is.null(hic$althalf)){

        if ((hic$half == "top" | hic$half == "bottom") && (half == "both")){

          stop("Invalid \'half\' of plot to annotate.", call. = FALSE)

        }

        if (hic$half == "top" & half == "bottom"){

          stop("Invalid \'half\' of plot to annotate.", call. = FALSE)

        }

        if (hic$half == "bottom" & half == "top"){

          stop("Invalid \'half\' of plot to annotate.", call. = FALSE)

        }

      } else {

        if (hic$althalf == "bottom"){

          if (!quiet) message(paste("Attempting to annotate pixels where",
                                    hic$chrom, "is on the x-axis and",
                                    hic$altchrom, "is on the y-axis."), call. = FALSE)

        } else if (hic$althalf == "top"){

          if (!quiet) message(paste("Attempting to annotate pixels where",
                                    hic$altchrom, "is on the x-axis and",
                                    hic$chrom, "is on the y-axis."), call. = FALSE)

        }

      }

    } else if (class(hic) == "bb_hicTriangle"
               | class(hic) == "bb_hicRectangle"){

      if (half == "both" | half == "bottom"){

        warning( paste0("Plot of class \'",
                        class(hic),
                        "\' detected. Pixels will automatically be annotated in the upper triangular of the plot."), call. = FALSE)

      }

    }

    ###### annotation #####

    ## Check type of annotation
    if (!type %in% c("box", "circle", "arrow")){

      stop("Invalid \'type\' of annotation.  Options are \'box\', \'circle\', or \'arrow\'.", call. = FALSE)

    }

  }

  ## Define a function that subsets loop data for hic region
  subset_loops <- function(hic, loops, object){

    ## chrom always in col1
    ## altchrom always in col4
    ## triangle hic plots will not have altchrom parameters
    if (class(hic) == "bb_hicTriangle" | class(hic) == "bb_hicRectangle"){
      loops_subset <- loops[which(loops[,1] == object$chrom
                                  & loops[,4] == object$chrom
                                  & loops[,2] >= object$chromstart
                                  & loops[,3] <= object$chromend
                                  & loops[,5] >= object$chromstart
                                  & loops[,6] <= object$chromend),]
    } else {
      loops_subset <- loops[which(loops[,1] == object$chrom
                                  & loops[,4] == object$altchrom
                                  & loops[,2] >= object$chromstart
                                  & loops[,3] <= object$chromend
                                  & loops[,5] >= object$altchromstart
                                  & loops[,6] <= object$altchromend),]
    }

    return(loops_subset)

  }

  ## Define a function that parses an inherited half
  inherit_half <- function(hic){

    if (class(hic) == "bb_hicSquare"){

      if (is.null(hic$althalf)){

        half <- hic$half

      } else {

        half <- hic$althalf

      }

    } else if (class(hic) == "bb_hicTriangle"
               | class(hic) == "bb_hicRectangle"){

      half <- "top"

    }

    return(half)

  }

  ## Define a function to add box annotation
  boxAnnotation <- function(df, hic, object, shift, half){

    side <- (as.numeric(df[6]) - as.numeric(df[5])) + (2 * shift * hic$resolution)

    if (half == "bottom"){

      center_x <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      rect1 <- rectGrob(x = center_x, y = center_y, width = side,
                        height = side, default.units = "native",
                        gp = object$gp)

      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = rect1), envir = bbEnv)

    } else if (half == "top"){

      center_x <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      rect1 <- rectGrob(x = center_x, y = center_y, width = side,
                        height = side, default.units = "native",
                        gp = object$gp)

      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = rect1), envir = bbEnv)

    } else if (half == "both"){

      ## BOTTOM
      center_x1 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y1 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))

      ## TOP
      center_x2 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y2 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))

      rect1 <- rectGrob(x = center_x1, y = center_y1, width = side,
                        height = side, default.units = "native",
                        gp = object$gp)
      rect2 <- rectGrob(x = center_x2, y = center_y2, width = side,
                        height = side, default.units = "native",
                        gp = object$gp)

      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = rect1), envir = bbEnv)
      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = rect2), envir = bbEnv)

    }

  }

  ## Define a function to add circle annotation
  circleAnnotation <- function(df, hic, object, shift, half){

    radius <- (0.5 * (as.numeric(df[6]) - as.numeric(df[5]))) + (shift * hic$resolution)

    if (half == "bottom"){

      center_x <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      circ1 <- circleGrob(x = center_x, y = center_y,
                          r = radius, default.units = "native",
                          gp = object$gp)

      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = circ1), envir = bbEnv)

    } else if (half == "top"){

      center_x <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      circ1 <- circleGrob(x = center_x, y = center_y,
                          r = radius, default.units = "native",
                          gp = object$gp)

      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = circ1), envir = bbEnv)

    } else if (half == "both"){

      ## BOTTOM
      center_x1 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y1 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))

      ## TOP
      center_x2 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y2 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))

      circ1 <- circleGrob(x = center_x1, y = center_y1,
                          r = radius, default.units = "native",
                          gp = object$gp)
      circ2 <- circleGrob(x = center_x2, y = center_y2,
                          r = radius, default.units = "native",
                          gp = object$gp)

      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = circ1), envir = bbEnv)
      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = circ2), envir = bbEnv)

    }

  }

  ## Define a function to add arrow annotation
  arrowAnnotation <- function(df, hic, object, shift, half){

    if (half == "bottom"){

     x0 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
     y0 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

     arrow1 <- segmentsGrob(x0 = x0, y0 = y0,
                            x1 = x0 + (shift * hic$resolution),
                            y1 = y0 - (shift * hic$resolution),
                            arrow = arrow(length = unit(0.1, "inches"),
                                          ends = "first",
                                          type = "closed"),
                            default.units = "native",
                            gp = object$gp)

     assign("loop_grobs",
            addGrob(gTree = get("loop_grobs", envir = bbEnv),
                    child = arrow1), envir = bbEnv)

    } else if (half == "top"){

      x0 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
      y0 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

      arrow1 <- segmentsGrob(x0 = x0, y0 = y0,
                             x1 = x0 - (shift * hic$resolution),
                             y1 = y0 + (shift * hic$resolution),
                             arrow = arrow(length = unit(0.1, "inches"),
                                           ends = "first",
                                           type = "closed"),
                             default.units = "native",
                             gp = object$gp)

      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = arrow1), envir = bbEnv)

    } else if (half == "both"){

      ## BOTTOM
      x01 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
      y01 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

      ## TOP
      x02 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
      y02 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

      arrow1 <- segmentsGrob(x0 = x01, y0 = y01,
                             x1 = x01 + (shift * hic$resolution),
                             y1 = y01 - (shift * hic$resolution),
                             arrow = arrow(length = unit(0.1, "inches"),
                                           ends = "first",
                                           type = "closed"),
                             default.units = "native",
                             gp = object$gp)
      arrow2 <- segmentsGrob(x0 = x02, y0 = y02,
                             x1 = x02 - (shift * hic$resolution),
                             y1 = y02 + (shift * hic$resolution),
                             arrow = arrow(length = unit(0.1, "inches"),
                                           ends = "first",
                                           type = "closed"),
                             default.units = "native",
                             gp = object$gp)

      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = arrow1), envir = bbEnv)
      assign("loop_grobs",
             addGrob(gTree = get("loop_grobs", envir = bbEnv),
                     child = arrow2), envir = bbEnv)

    }

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(half)) half <- NULL
  if(missing(shift)) shift <- NULL
  if(missing(type)) type <- NULL
  if(missing(quiet)) quiet <- NULL

  ## Check if hic/loops arguments are missing (could be in object)
  if(!hasArg(plot)) plot <- NULL
  if(!hasArg(data)) data <- NULL

  ## Compile all parameters into an internal object
  bb_loopsInternal <- structure(list(plot = plot, data = data, half = half,
                                     shift = shift, type = type, quiet = quiet,
                                     gp = gpar()), class = "bb_loopsInternal")

  bb_loopsInternal <- parseParams(bb_params = params,
                                  object_params = bb_loopsInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_loopsInternal$half)) bb_loopsInternal$half <- "inherit"
  if(is.null(bb_loopsInternal$shift)) bb_loopsInternal$shift <- 4
  if(is.null(bb_loopsInternal$type)) bb_loopsInternal$type <- "box"
  if(is.null(bb_loopsInternal$quiet)) bb_loopsInternal$quiet <- FALSE

  ## Set gp
  bb_loopsInternal$gp <- setGP(gpList = bb_loopsInternal$gp,
                               params = bb_loopsInternal, ...)

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT: GET REGION/DIMENSIONS FROM HIC PLOT INPUT
  # ======================================================================================================================================================================================

  bb_loops <- structure(list(chrom = bb_loopsInternal$plot$chrom,
                             chromstart = bb_loopsInternal$plot$chromstart,
                             chromend = bb_loopsInternal$plot$chromend,
                             altchrom = bb_loopsInternal$plot$altchrom,
                             altchromstart = bb_loopsInternal$plot$altchromstart,
                             altchromend = bb_loopsInternal$plot$altchromend,
                             assembly = bb_loopsInternal$plot$assembly,
                             x = bb_loopsInternal$plot$x,
                             y = bb_loopsInternal$plot$y,
                             width = bb_loopsInternal$plot$width,
                             height = bb_loopsInternal$plot$height,
                             just = bb_loopsInternal$plot$just, grobs = NULL),
                        class = "bb_pixel")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot annotate Hi-C pixels without a BentoBox page.")
  if(is.null(bb_loopsInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_loopsInternal$data)) stop("argument \"data\" is missing, with no default.", call. = FALSE)


  errorcheck_bb_annoLoops(hic = bb_loopsInternal$plot,
                          loops = bb_loopsInternal$data,
                          half = bb_loopsInternal$half,
                          type = bb_loopsInternal$type,
                          quiet = bb_loopsInternal$quiet)

  # ======================================================================================================================================================================================
  # PARSE INHERITED HALF
  # ======================================================================================================================================================================================

  half <- bb_loopsInternal$half
  if (half == "inherit"){

    half <- inherit_half(hic = bb_loopsInternal$plot)

  }

  if (class(bb_loopsInternal$plot) == "bb_hicTriangle"
      | class(bb_loopsInternal$plot) == "bb_hicRectangle"){

    half <- "top"

  }

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  loops <- bb_loopsInternal$data
  if (!"data.frame" %in% class(loops)){

    loops <- as.data.frame(data.table::fread(loops))
    if (nrow(loops) < 1){
      warning("\'data\' input contains no values.", call. = FALSE)
    }


  }

  ## Check format of chromosomes in columns 1 and 4
  if (bb_loops$assembly$Genome == "hg19"){

    checkChr <- function(chr){
      return(grepl("chr", chr))
    }

    col1Checks <- unlist(lapply(loops[,1], checkChr))
    if (any(col1Checks == FALSE)){
      stop("Chromosomes in column 1 are in invalid format for hg19 genome assembly. Please specify chromosomes as a string with the following format: 'chr1'.", call. = FALSE)
    }
    col4Checks <- unlist(lapply(loops[,4], checkChr))
    if (any(col4Checks == FALSE)){
      stop("Chromosomes in column 4 are in invalid format for hg19 genome assembly. Please specify chromosomes as a string with the following format: 'chr1'.", call. = FALSE)
    }


  }

  # ======================================================================================================================================================================================
  # SUBSET FOR LOOPS IN REGION
  # ======================================================================================================================================================================================
  ## Assuming loops are in first six columns only
  loops <- loops[,c(seq(1,6))]

  loops_subset <- subset_loops(hic = bb_loopsInternal$plot, loops = loops,
                               object = bb_loops)
  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_pixel",
                    length(grep(pattern = "bb_pixel",
                                x = currentViewports)) + 1)

  ## Make viewport based on hic input viewport
  if (class(bb_loopsInternal$plot) == "bb_hicSquare"){

    vp <- viewport(height = bb_loopsInternal$plot$grobs$vp$height,
                   width = bb_loopsInternal$plot$grobs$vp$width,
                   x = bb_loopsInternal$plot$grobs$vp$x,
                   y = bb_loopsInternal$plot$grobs$vp$y,
                   clip = "on",
                   xscale = bb_loopsInternal$plot$grobs$vp$xscale,
                   yscale = bb_loopsInternal$plot$grobs$vp$yscale,
                   just = bb_loopsInternal$plot$grobs$vp$justification,
                   name = vp_name)

  } else if (class(bb_loopsInternal$plot) == "bb_hicTriangle"){


    width <- convertUnit(bb_loopsInternal$plot$outsideVP$width,
                         unitTo = get("page_units", bbEnv), valueOnly = TRUE)

    vp <- viewport(height = unit(width/sqrt(two), get("page_units", bbEnv)),
                   width = unit(width/sqrt(two), get("page_units", bbEnv)),
                   x = bb_loopsInternal$plot$outsideVP$x,
                   y = bb_loopsInternal$plot$outsideVP$y,
                   xscale = bb_loopsInternal$plot$grobs$vp$xscale,
                   yscale = bb_loopsInternal$plot$grobs$vp$yscale,
                   just = bb_loopsInternal$plot$outsideVP$justification,
                   name = vp_name,
                   angle = -45)

  } else if (class(bb_loopsInternal$plot) == "bb_hicRectangle"){

    side <- convertUnit(bb_loopsInternal$plot$grobs$vp$width,
                        unitTo = get("page_units", bbEnv))

    ## Get bottom left coord of outsideVP
    bottomLeft <- vp_bottomLeft(viewport = bb_loopsInternal$plot$outsideVP)

    ## Convert adjusted chromstart to page units within outsideVP and add to bottomLeft x
    seekViewport(name = bb_loopsInternal$plot$outsideVP$name)
    xCoord <- convertX(unit(bb_loopsInternal$plot$grobs$vp$x),
                       unitTo = get("page_units", bbEnv)) + bottomLeft[[1]]
    seekViewport(name = "bb_page")

    vp <- viewport(height = side, width = side,
                   x = xCoord, y = bottomLeft[[2]],
                   xscale = bb_loopsInternal$plot$grobs$vp$xscale,
                   yscale = bb_loopsInternal$plot$grobs$vp$yscale,
                   just = c("left", "bottom"),
                   name = vp_name,
                   angle = -45)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE OF GROBS
  # ======================================================================================================================================================================================

  assign("loop_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  if (nrow(loops_subset) > 0){

    if (bb_loopsInternal$type == "box"){

      bb_loopsInternal$gp$fill <- NA
      invisible(apply(loops_subset, 1, boxAnnotation,
                      hic = bb_loopsInternal$plot, object = bb_loopsInternal,
                      shift = bb_loopsInternal$shift, half = half))

    } else if (bb_loopsInternal$type == "circle"){
      bb_loopsInternal$gp$fill <- NA
      invisible(apply(loops_subset, 1, circleAnnotation,
                      hic = bb_loopsInternal$plot, object = bb_loopsInternal,
                      shift = bb_loopsInternal$shift, half = half))

    } else if (bb_loopsInternal$type == "arrow"){
      if (is.null(bb_loopsInternal$gp$col)
          & is.null(bb_loopsInternal$gp$fill)){
        bb_loopsInternal$gp$fill <- "black"
      } else {
        if(is.null(bb_loopsInternal$gp$fill)){
          bb_loopsInternal$gp$fill <- bb_loopsInternal$gp$col
        }
      }

      invisible(apply(loops_subset, 1, arrowAnnotation,
                      hic = bb_loopsInternal$plot, object = bb_loopsInternal,
                      shift = bb_loopsInternal$shift, half = half))

    }

  } else {

    warning("No pixels found in region.", call. = FALSE)
  }


  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  bb_loops$grobs <- get("loop_grobs", envir = bbEnv)
  grid.draw(bb_loops$grobs)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_pixel[", vp_name, "]"))
  invisible(bb_loops)

}
