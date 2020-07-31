#' plots data stored in BEDPE format in a "ribbon arch" style
#'
#' @param bedpe bed paired end data to be plotted
#' @param chrom chromsome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param loopHeight single value or vector specifying the arch heights;
#' when NULL, all arches will be the same height, filling up the given viewport
#' @param style style of arches: "3D" or "2D"
#' @param arch arch curvature, determined by number of points along the curve
#' @param position position of arches above or below the x-axis: "top" or "bottom"
#' @param fillcolor single value, vector, or palette specifying colors of arches
#' @param linecolor linecolor of line outlining arches
#' @param colorby name of column in bedpe data to scale colors by
#' @param colorbyrange the range of values to apply a colorby palette scale to, if colorby values are numeric
#' @param assembly desired genome assembly
#' @param transp transparency of arches
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param just A string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @return Function will plot a ribbon arches plot and return a bb_arches object
#'
#' @export

bb_plotRibbonArches <- function(bedpe, chrom, chromstart = NULL, chromend = NULL, loopHeight = NULL, style = "3D", arch = 5, position = "top", fillcolor = "lightgrey",
                                linecolor = NA, colorby = NULL, colorbyrange = NULL, assembly = "hg19", transp = 0.4, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"),
                                default.units = "inches", draw = T){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors
  errorcheck_bbArches <- function(bedpe, arches_plot, style, position, colorby){
    ## Can't have only one NULL chromstart or chromend
    if ((is.null(arches_plot$chromstart) & !is.null(arches_plot$chromend)) | (is.null(arches_plot$chromend) & !is.null(arches_plot$chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.", call. = FALSE)

    }


    if (!is.null(arches_plot$chromstart) & !is.null(arches_plot$chromend)){

      ## chromend > chromstart
      if (arches_plot$chromend < arches_plot$chromstart){

        stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)


      }

    }

    if (!style %in% c("3D", "2D")){
      stop("Invalid \'style\' input. Options are \'3D\' and \'2D\'.", call. = FALSE)
    }

    if (!position %in% c("top", "bottom")){
      stop("Invalid \'position\' input. Options are \'top\' and \'bottom\'.", call. = FALSE)
    }

    if (!is.null(colorby)){
      if (!any(colnames(bedpe) == colorby)){
        stop("Colorby column not found in data. Check colorby column name.", call. = FALSE)
      }

      if (length(which(colnames(bedpe) == colorby)) > 1){
        stop("Multiple matching colorby columns found in data. Please provide colorby column name with only one occurrence.", call. = FALSE)
      }
    }
}

  ## Define a function that normalizes arch heights
  normHeights <- function(height, position, min, max){

    ## First normalize from 0 to 1
    newHeight <- (height - min)/(max-min)

    ## Then scale to [0, 1.38] for top position and [-0.38, 1] for bottom position
    if (position == "top"){
      finalHeight <- newHeight * 1.38
    } else {
      finalHeight <- (newHeight * 1.38) + (-0.38)
    }

    return(finalHeight)

  }

  ## Define a function that creates ribbon arch grobs
  drawRibbons <- function(df, style, arch, position, linecolor, transp){

    x1 <- df[[2]]
    x2 <- df[[3]]
    y1 <- df[[5]]
    y2 <- df[[6]]
    fillCol <- df[[7]]
    outerHeight <- as.numeric(df[[9]])


    if (style == "3D"){
      x1 <- df[[3]]
      x2 <- df[[2]]
    }

    ## Designate bezier control points
    innerX = unit(seq(x2, y1, length.out = arch)[c(1:2, (arch-1):arch)], "native")
    outerX = unit(seq(x1, y2, length.out = arch)[c(1:2, (arch-1):arch)], "native")

    if (position == "top"){
      innerHeight <- outerHeight - 0.01
      ## Switch y-positions for top plotting
      innerY = unit(c(0, innerHeight, innerHeight, 0), "npc")
      outerY = unit(c(0, outerHeight, outerHeight, 0), "npc")

    } else {
      innerHeight <- outerHeight + 0.01
      ## Switch y-positions for bottom plotting
      innerY = unit(c(1, innerHeight, innerHeight, 1), "npc")
      outerY = unit(c(1, outerHeight, outerHeight, 1), "npc")
    }

    ## Calculate loop arcs using bezier curves
    innerLoop <- bezierGrob(x = innerX, y = innerY)
    outerLoop <- bezierGrob(x = outerX, y = outerY)

    ## Extract points from bezier curves
    innerBP <- bezierPoints(innerLoop)
    outerBP <- bezierPoints(outerLoop)

    ## Connect points, convert to proper units and draw polygons
    archGrob <- polygonGrob(x = unit(c(convertX(outerBP$x, "native"), rev(convertX(innerBP$x, "native"))), "native"),
                 y = unit(c(convertY(outerBP$y, "npc"), rev(convertY(innerBP$y, "npc"))), "npc"),
                 gp = gpar(fill = fillCol, col = linecolor, alpha = transp))
    assign("arches_grobs", addGrob(gTree = get("arches_grobs", envir = bbEnv), child = archGrob), envir = bbEnv)


  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  arches_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, zrange = NULL, color_palette = NULL, width = width,
                                height = height, x = x, y = y, justification = just, grobs = NULL, assembly = assembly), class = "bb_pileup")
  attr(x = arches_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CHECK PLACEMENT
  # ======================================================================================================================================================================================

  check_placement(object = arches_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  arches_plot <- defaultUnits(object = arches_plot, default.units = default.units)

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if ("data.frame" %in% class(bedpe)){
    bedpe <- as.data.frame(bedpe)
  } else {
    bedpe <- as.data.frame(data.table::fread(bedpe))
  }

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bbArches(bedpe = bedpe, arches_plot = arches_plot, style = style, position = position, colorby = colorby)

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  if (is.null(chromstart) & is.null(chromend)){

    if (assembly == "hg19"){
      genome <- bb_hg19
    }

    arches_plot$chromstart <- 1
    arches_plot$chromend <- genome[which(genome$chrom == chrom),]$length

  }

  bedpe <- bedpe[which(bedpe[,1] == arches_plot$chrom & bedpe[,4] == arches_plot$chrom
                       & bedpe[,2] >= arches_plot$chromstart & bedpe[,3] <= arches_plot$chromend
                              & bedpe[,5] >= arches_plot$chromstart & bedpe[,6] <= arches_plot$chromend),]

  # ======================================================================================================================================================================================
  # COLORBY AND COLORS
  # ======================================================================================================================================================================================

   if (!is.null(colorby)){
    colorbyCol <- which(colnames(bedpe) == colorby)
    colorbyCol <- bedpe[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      colorbyCol <- as.numeric(colorbyCol)
    } else {

      if (is.null(colorbyrange)){
        colorbyrange <- c(min(colorbyCol), max(colorbyCol))
        arches_plot$zrange <- colorbyrange
      }

    }

    if (class(fillcolor) == "function"){
      colors <- bb_maptocolors(colorbyCol, fillcolor, range = colorbyrange)
      sorted_colors <- unique(colors[order(colorbyCol)])
      arches_plot$color_palette <- sorted_colors

    } else {
      colorbyColfac <- factor(colorbyCol)
      mappedColors <- rep(fillcolor, ceiling(length(levels(colorbyColfac))/length(fillcolor)))
      colors <- mappedColors[colorbyCol]

    }

  } else {
    if (class(fillcolor) == "function"){
      colors <- fillcolor(nrow(bedpe))
    } else {
      if (length(fillcolor) == 1){
        colors <- rep(fillcolor, nrow(bedpe))
      } else {
        colors <- rep(fillcolor, ceiling(nrow(bedpe)/length(fillcolor)))[1:nrow(bedpe)]
      }
    }

  }

  bedpe <- bedpe[,c(1:6)]
  bedpe$color <- colors

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_arches", length(grep(pattern = "bb_arches", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(0.5, "npc"), width = unit(1, "npc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   xscale = c(arches_plot$chromstart, arches_plot$chromend),
                   clip = "on",
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      vp$name <- "bb_arches1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = arches_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   xscale = c(arches_plot$chromstart, arches_plot$chromend),
                   clip = "on",
                   just = just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # HEIGHTS
  # ======================================================================================================================================================================================

  if (is.null(loopHeight)){
    bedpe$height <- rep(1, nrow(bedpe))
  } else if (length(loopHeight) == 1){
    bedpe$height <- rep(loopHeight, nrow(bedpe))
  } else {
    bedpe$height <- loopHeight[1:nrow(bedpe)]
  }

  if (length(bedpe$height) > 0){
    bedpe$normHeight <- lapply(bedpe$height, normHeights, position = position, min = 0, max = max(bedpe$height))
  }


  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("arches_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # GROBS
  # ======================================================================================================================================================================================

  if (nrow(bedpe) > 0){

    invisible(apply(bedpe, 1, drawRibbons, style = style, arch = arch, position = position,
                    linecolor = linecolor, transp = transp))

  } else {
    warning("Bedpe contains no values.", call. = FALSE)
  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("arches_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  arches_plot$grobs <-  get("arches_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(arches_plot)

}
