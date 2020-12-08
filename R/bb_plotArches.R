#' plots data stored in BEDPE format in a "ribbon arch" style
#'
#' @param bedpe bed paired end data to be plotted
#' @param chrom chromsome of region to be plotted
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param chromstart start position
#' @param chromend end position
#' @param cutoffArches a logical value indicating whether to include arches that get cutoff in the region
#' @param archHeight single value or vector specifying the arch heights;
#' when NULL, all arches will be the same height, filling up the given viewport
#' @param style style of arches: "2D" or "3D"
#' @param arch arch curvature, determined by number of points along the curve
#' @param position position of arches above or below the x-axis: "top" or "bottom"
#' @param fillcolor single value, vector, or palette specifying colors of arches
#' @param linecolor linecolor of line outlining arches
#' @param colorby name of column in bedpe data to scale colors by
#' @param colorbyrange the range of values to apply a colorby palette scale to, if colorby values are numeric
#' @param assembly desired genome assembly
#' @param alpha transparency of arches
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

bb_plotArches <- function(bedpe, chrom, params = NULL, chromstart = NULL, chromend = NULL, cutoffArches = FALSE, archHeight = NULL, style = "2D", arch = 5, position = "top",
                                fillcolor = "lightgrey", linecolor = NA, colorby = NULL, colorbyrange = NULL, assembly = "hg19", alpha = 0.4, x = NULL, y = NULL, width = NULL,
                                height = NULL, just = c("left", "top"), default.units = "inches", draw = T, ...){

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

  ## Define a function that will produce a yscale for arches based on height/position
  height_yscale <- function(heights, position){

    if (length(heights) == 1){

      if (position == "top"){
        yscale <- c(0, heights)
      } else {
        yscale <- c(heights, 0)
      }

    } else {

      if (position == "top"){
        yscale <- c(0, max(heights))
      } else {
        yscale <- c(max(heights), 0)
      }


    }
    return(yscale)
  }

  ## Define a function that normalizes arch heights
  normHeights <- function(height, position, min, max){

    ## First normalize from 0 to 1
    newHeight <- (height - min)/(max-min)

    ## Then scale to a range of 1.38
    finalHeight <- newHeight * 1.38

    return(finalHeight)

  }

  ## Define a function that creates ribbon arch grobs
  drawRibbons <- function(df, style, arch, position, linecolor, transp, ...){

    x1 <- df[[2]]
    x2 <- df[[3]]
    y1 <- df[[5]]
    y2 <- df[[6]]
    fillCol <- df[[7]]
    outerHeight <- as.numeric(df[[9]])
    innerHeight <- outerHeight - 0.01

    if (style == "3D"){
      x1 <- df[[3]]
      x2 <- df[[2]]
    }

    ## Designate bezier control points
    innerX = unit(seq(x2, y1, length.out = arch)[c(1:2, (arch-1):arch)], "native")
    outerX = unit(seq(x1, y2, length.out = arch)[c(1:2, (arch-1):arch)], "native")

    if (position == "top"){
      ## Switch y-positions for top plotting
      innerY = unit(c(0, innerHeight, innerHeight, 0), "npc")
      outerY = unit(c(0, outerHeight, outerHeight, 0), "npc")

    } else {
      ## Switch y-positions for bottom plotting
      innerY = unit(c(1, 1-innerHeight, 1-innerHeight, 1), "npc")
      outerY = unit(c(1, 1-outerHeight, 1-outerHeight, 1), "npc")
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
                 gp = gpar(fill = fillCol, col = linecolor, alpha = transp, ...))
    assign("arches_grobs", addGrob(gTree = get("arches_grobs", envir = bbEnv), child = archGrob), envir = bbEnv)


  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(cutoffArches)) cutoffArches <- NULL
  if(missing(style)) style <- NULL
  if(missing(arch)) arch <- NULL
  if(missing(position)) position <- NULL
  if(missing(fillcolor)) fillcolor <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(alpha)) alpha <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if bedpe/chrom arguments are missing (could be in object)
  if(!hasArg(bedpe)) bedpe <- NULL
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_archInternal <- structure(list(bedpe = bedpe, chrom = chrom, chromstart = chromstart, chromend = chromend, cutoffArches = cutoffArches, archHeight = archHeight, style = style,
                                    arch = arch, position = position, fillcolor = fillcolor, linecolor = linecolor, colorby = colorby, colorbyrange = colorbyrange,
                                    assembly = assembly, alpha = alpha, x = x, y = y, width = width, height = height, just = just,
                                    default.units = default.units, draw = draw), class = "bb_archInternal")

  bb_archInternal <- parseParams(bb_params = params, object_params = bb_archInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_archInternal$cutoffArches)) bb_archInternal$cutoffArches <- FALSE
  if(is.null(bb_archInternal$style)) bb_archInternal$style <- "2D"
  if(is.null(bb_archInternal$arch)) bb_archInternal$arch <- 5
  if(is.null(bb_archInternal$position)) bb_archInternal$position <- "top"
  if(is.null(bb_archInternal$fillcolor)) bb_archInternal$fillcolor <- "lightgrey"
  if(is.null(bb_archInternal$linecolor)) bb_archInternal$linecolor <- NA
  if(is.null(bb_archInternal$assembly)) bb_archInternal$assembly <- "hg19"
  if(is.null(bb_archInternal$alpha)) bb_archInternal$alpha <- 0.4
  if(is.null(bb_archInternal$just)) bb_archInternal$just <- c("left", "top")
  if(is.null(bb_archInternal$default.units)) bb_archInternal$default.units <- "inches"
  if(is.null(bb_archInternal$draw)) bb_archInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  arches_plot <- structure(list(chrom = bb_archInternal$chrom, chromstart = bb_archInternal$chromstart, chromend = bb_archInternal$chromend, zrange = NULL, color_palette = NULL,
                                width = bb_archInternal$width, height = bb_archInternal$height, x = bb_archInternal$x, y = bb_archInternal$y, justification = bb_archInternal$just,
                                grobs = NULL, assembly = bb_archInternal$assembly), class = "bb_arches")
  attr(x = arches_plot, which = "plotted") <- bb_archInternal$draw

  # ======================================================================================================================================================================================
  # CHECK PLACEMENT/ARGUMENT ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_archInternal$bedpe)) stop("argument \"bedpe\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_archInternal$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
  check_placement(object = arches_plot)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  arches_plot$assembly <- parse_bbAssembly(assembly = arches_plot$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  arches_plot <- defaultUnits(object = arches_plot, default.units = bb_archInternal$default.units)

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if ("data.frame" %in% class(bb_archInternal$bedpe)){
    bedpe <- as.data.frame(bb_archInternal$bedpe)
  } else {
    bedpe <- as.data.frame(data.table::fread(bb_archInternal$bedpe))
  }

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bbArches(bedpe = bedpe, arches_plot = arches_plot, style = bb_archInternal$style, position = bb_archInternal$position, colorby = bb_archInternal$colorby)

  # ======================================================================================================================================================================================
  # WHOLE CHROM DATA AND XSCALE
  # ======================================================================================================================================================================================

  if (is.null(arches_plot$chromstart) & is.null(arches_plot$chromend)){

    txdbChecks <- check_loadedPackage(package = arches_plot$assembly$TxDb, message = paste(paste0("`", arches_plot$assembly$TxDb,"`"),
                                                                                          "not loaded. Please install and load to plot full chromosome ribbon arches."))
    xscale <- c(0, 1)
    if (txdbChecks == TRUE){
      tx_db <- eval(parse(text = arches_plot$assembly$TxDb))
      assembly_data <- seqlengths(tx_db)

      if (!arches_plot$chrom %in% names(assembly_data)){
        txdbChecks <- FALSE
        warning(paste("Chromosome", paste0("'", arches_plot$chrom, "'"), "not found in", paste0("`", arches_plot$assembly$TxDb, "`"), "and data for entire chromosome cannot be plotted."), call. = FALSE)
      } else {
        arches_plot$chromstart <- 1
        arches_plot$chromend <- assembly_data[[arches_plot$chrom]]
        xscale <- c(arches_plot$chromstart, arches_plot$chromend)

      }

    }

  } else {
    txdbChecks <- TRUE
    xscale <- c(arches_plot$chromstart, arches_plot$chromend)
  }

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  if (!is.null(arches_plot$chromstart) & !is.null(arches_plot$chromend)){

    if (bb_archInternal$cutoffArches == FALSE){
      bedpe <- bedpe[which(bedpe[,1] == arches_plot$chrom & bedpe[,4] == arches_plot$chrom
                           & bedpe[,2] >= arches_plot$chromstart & bedpe[,3] <= arches_plot$chromend
                           & bedpe[,5] >= arches_plot$chromstart & bedpe[,6] <= arches_plot$chromend),]
    } else {
      bedpe <- bedpe[bedpe[[1]] == arches_plot$chrom & bedpe[[4]] == arches_plot$chrom &
                       ((bedpe[[3]] >= arches_plot$chromstart & bedpe[[3]] <= arches_plot$chromend) |
                          (bedpe[[5]] <= arches_plot$chromstart & bedpe[[5]] >= arches_plot$chromend))]
    }

  } else {

    bedpe <- data.frame(matrix(nrow = 0, ncol = 6))
  }

  # ======================================================================================================================================================================================
  # COLORBY AND COLORS
  # ======================================================================================================================================================================================

   if (!is.null(bb_archInternal$colorby) & nrow(bedpe) > 0){
    colorbyCol <- which(colnames(bedpe) == bb_archInternal$colorby)
    colorbyCol <- bedpe[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      colorbyCol <- as.numeric(colorbyCol)
    }

    if (is.null(bb_archInternal$colorbyrange)){
      colorbyrange <- c(min(colorbyCol), max(colorbyCol))
      arches_plot$zrange <- colorbyrange
    }

    if (class(bb_archInternal$fillcolor) == "function"){
      colors <- bb_maptocolors(colorbyCol, bb_archInternal$fillcolor, range = arches_plot$zrange)
      sorted_colors <- unique(colors[order(colorbyCol)])
      arches_plot$color_palette <- sorted_colors

    } else {
      colorbyColfac <- factor(colorbyCol)
      mappedColors <- rep(bb_archInternal$fillcolor, ceiling(length(levels(colorbyColfac))/length(bb_archInternal$fillcolor)))
      colors <- mappedColors[colorbyCol]

    }

  } else {
    if (class(bb_archInternal$fillcolor) == "function"){
      colors <- bb_archInternal$fillcolor(nrow(bedpe))
    } else {
      if (length(bb_archInternal$fillcolor) == 1){
        colors <- rep(bb_archInternal$fillcolor, nrow(bedpe))
      } else {
        colors <- rep(bb_archInternal$fillcolor, ceiling(nrow(bedpe)/length(bb_archInternal$fillcolor)))[1:nrow(bedpe)]
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
  if (is.null(arches_plot$x) & is.null(arches_plot$y)){

    vp <- viewport(height = unit(0.5, "npc"), width = unit(1, "npc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   xscale = xscale,
                   clip = "on",
                   just = "center",
                   name = vp_name)

    if (bb_archInternal$draw == TRUE){

      vp$name <- "bb_arches1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = arches_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   xscale = xscale,
                   clip = "on",
                   just = bb_archInternal$just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # HEIGHTS
  # ======================================================================================================================================================================================

  if (is.null(bb_archInternal$archHeight)){
    bedpe$height <- rep(1, nrow(bedpe))
  } else if (length(bb_archInternal$archHeight) == 1){
    bedpe$height <- rep(bb_archInternal$archHeight, nrow(bedpe))
    yscale <- height_yscale(heights = bb_archInternal$archHeight, position = bb_archInternal$position)
    vp$yscale <- yscale
  } else {
    bedpe$height <- bb_archInternal$archHeight[1:nrow(bedpe)]
    yscale <- height_yscale(heights = bb_archInternal$archHeight, position = bb_archInternal$position)
    vp$yscale <- yscale
  }

  if (length(bedpe$height) > 0){
    bedpe$normHeight <- lapply(bedpe$height, normHeights, position = bb_archInternal$position, min = 0, max = max(bedpe$height))
  }


  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("arches_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # GROBS
  # ======================================================================================================================================================================================

  if (nrow(bedpe) > 0){

    invisible(apply(bedpe, 1, drawRibbons, style = bb_archInternal$style, arch = bb_archInternal$arch, position = bb_archInternal$position,
                    linecolor = bb_archInternal$linecolor, transp = bb_archInternal$alpha, ...))

  } else {

    if (txdbChecks == TRUE){
      warning("Bedpe contains no values.", call. = FALSE)
    }

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_archInternal$draw == TRUE){

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
