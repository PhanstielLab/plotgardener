#' plots paired-end data
#'
#' @param bedpe bed paired end data to be plotted
#' @param chrom chromsome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param fillcolor single value or vector specifying colors of bedpe elements
#' @param colorby vector to scale colors by
#' @param colorbycol palette to apply color scale to (only valid when colorby is not NULL)
#' @param colorbyrange the range of values to apply the color scale to
#' @param linecolor border color
#' @param boxHeight height of boxes at either end of bedpe element
#' @param spaceHeight height of space between boxes of different bedpe elements
#' @param limitDots logical value indicating whether to plot "..." to indicate additional, unplotted bedpe elements
#' @param x A unit object specifying x-location
#' @param y A unit object specifying y-location
#' @param width A unit object specifying width
#' @param height A unit object specifying height
#' @param just string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @export

bb_plotBedpe <- function(bedpe, chrom, chromstart, chromend, fillcolor = "black", colorby = NULL, colorbycol = NULL, colorbyrange = NULL, linecolor = NULL, boxHeight = unit(0.025, "inches"),
                         spaceHeight = unit(.025, "inches"), limitDots = TRUE, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotBedpe <- function(bedpe, bedpe_plot){



    ## bedpe file
    ## if it's a file path, it needs to exist
    if (!"data.frame" %in% class(bedpe)){

      # ## File extension
      if (file_ext(loops) != "bedpe"){

        stop("Invalid input. File must have a \".bedpe\" extension")

      }

      ## File existence
      if (!file.exists(loops)){

        stop(paste("File", loops, "does not exist."), call. = FALSE)

      }

    }


    ## chrom needs to be a character
    if (class(bedpe_plot$chrom) != "character"){

      stop("Invalid \'chrom\'; input must be a string of form \"chr\"_.", call. = FALSE)

    }

    ## chromend > chromstart
    if (bedpe_plot$chromend < bedpe_plot$chromstart){

      stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)


    }

  }


  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bedpe_plot <- structure(list(chrom = chrom, chromstart = as.numeric(chromstart), chromend = as.numeric(chromend), fillcolor = fillcolor,
                                 linecolor = linecolor, colorby = colorby, colorbycol = colorbycol, colorbyrange = colorbyrange,
                                 width = width, height = height, x = x, y = y, justification = just, grobs = NULL), class = "bb_bedpe")
  attr(x = bedpe_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = bedpe_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  if (!(is.null(x) & is.null(y))){

    if (is.numeric(class(x)) | is.numeric(class(y)) | is.numeric(class(width)) | is.numeric(class(height))){

      if (is.null(default.units)){

        stop("One or more placement coordinates detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      bedpe_plot$x <- unit(x, default.units)
      bedpe_plot$y <- unit(y, default.units)
      bedpe_plot$width <- unit(width, default.units)
      bedpe_plot$height <- unit(height, default.units)

    }

  }

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if (!"data.frame" %in% class(bedpe)){

    bedpe <- data.table::fread(bedpe)

  }

  # ======================================================================================================================================================================================
  # ORGANIZE DATA
  # ======================================================================================================================================================================================

  ## Assuming data is in first six columns only
  bedpe <- bedpe[,1:6]
  colnames(bedpe) = c("chrom1", "start1", "stop1", "chrom2", "start2", "stop2")

  ## Get appropriate starts/stops
  start1 <- apply(bedpe[,c(2,3)], 1, min)
  stop1 <-  apply(bedpe[,c(2,3)], 1, max)
  start2 <- apply(bedpe[,c(5,6)], 1, min)
  stop2 <-  apply(bedpe[,c(5,6)], 1, max)
  bedpe$start1 <- start1
  bedpe$stop1 <- stop1
  bedpe$start2 <- start2
  bedpe$stop2 <- stop2

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  bedpe <- bedpe[which(bedpe[,1] == chrom & bedpe[,4] == chrom & bedpe[,2] >= chromstart & bedpe[,3] <= chromend
                       & bedpe[,5] >= chromstart & bedpe[,6] <= chromend),]

  # ======================================================================================================================================================================================
  # GET BOX WIDTHS
  # ======================================================================================================================================================================================
  bedpe$width1 <- bedpe$stop1 - bedpe$start1
  bedpe$width2 <- bedpe$stop2 - bedpe$start2
  # ======================================================================================================================================================================================
  # COLORS
  # ======================================================================================================================================================================================

  if (length(fillcolor) == 1){

    bedpe$color <- rep(fillcolor, nrow(bedpe))

  } else if (length(fillcolor) > 1){

    bedpe$color <- fillcolor
  }

  ## Add colorby column
  if (!is.null(colorby)){

    bedpedata$colorbyvalue <- colorby
    bedpe$color <- bb_maptocolors(bedpe$colorbyvalue, colorbycol, range = colorbyrange)

    if (is.null(colorbyrange)){

      colorbyrange <- c(min(bedpe$colorbyvalue), max(bedpe$colorbyvalue))

    }
  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_bedpe", length(grep(pattern = "bb_bedpe", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(1, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = c(chromstart, chromend),
                   yscale = c(0, 1),
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      vp$name <- "bb_bedpe1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = bedpe_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(chromstart, chromend),
                   yscale = c(0, convertHeight(height, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)),
                   just = just,
                   name = vp_name)
  }



  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("bedpe_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (nrow(bedpe) > 0){

    ## Determine how many bepe elements are going to fit based on boxHeight and space
    if (is.null(x) & is.null(y)){

      pushViewport(vp)
      boxHeight <- convertHeight(boxHeight, unitTo = "npc", valueOnly = T)
      spaceHeight <- convertHeight(spaceHeight, unitTo = "npc", valueOnly = T)
      upViewport()

    } else {

      boxHeight <- convertHeight(boxHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
      spaceHeight <- convertHeight(spaceHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

    }

    limit <- floor((vp$yscale[2] + spaceHeight)/(boxHeight + spaceHeight))

    if (nrow(bedpe) > limit){
      bedpe <- bedpe[1:limit,]
      warning("Not enough plotting space for all provided bedpe elements.", call. = FALSE)
    } else if (nrow(bedpe) < limit){
      warning("Excess plotting space for the provided bedpe elements.", call. = FALSE)
    }

    bedpe$rowNum <- seq(0, nrow(bedpe) - 1, 1)
    bedpe$y <- vp$yscale[2] - (bedpe$rowNum*(boxHeight + spaceHeight))

    bedpeRect1 <- rectGrob(x = bedpe$start1,
                           y = bedpe$y,
                           width = bedpe$width1,
                           height = boxHeight,
                           just = c("left", "top"),
                           default.units = "native",
                           gp = gpar(fill = bedpe$color, col = linecolor))

    bedpeRect2 <- rectGrob(x = bedpe$start2,
                           y = bedpe$y,
                           width = bedpe$width2,
                           height = boxHeight,
                           just = c("left", "top"),
                           default.units = "native",
                           gp = gpar(fill = bedpe$color, col = linecolor))

    bedpeLine <- segmentsGrob(x0 = bedpe$stop1,
                              y0 = bedpe$y - 0.5*boxHeight,
                              x1 = bedpe$start2,
                              y1 = bedpe$y - 0.5*boxHeight,
                              default.units = "native",
                              gp = gpar(col = bedpe$color))


    assign("bedpe_grobs", setChildren(get("bedpe_grobs", envir = bbEnv), children = gList(bedpeRect1, bedpeRect2, bedpeLine)), envir = bbEnv)

    if (limitDots == TRUE){

      dotsGrob <- textGrob(label = "...", x = unit(1, "npc"), y = unit(0, "npc"), just = c("right", "bottom"))

      assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = dotsGrob), envir = bbEnv)

    }

  } else {

    ## Just make a rectangle
    bedpeGrob <- rectGrob()
    assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = bedpeGrob), envir = bbEnv)
    warning("Bedpe contains no values.", call. = FALSE)

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("bedpe_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  bedpe_plot$grobs <-  get("bedpe_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bedpe_plot)

}
