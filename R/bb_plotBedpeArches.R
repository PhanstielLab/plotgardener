#' Plot paired-end BEDPE data in an arch style
#'
#' @param data A string specifying the BEDPE file path or a dataframe
#' in BEDPE format specifying data to be plotted.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[BentoBox]{bb_assembly} object.
#' Default value is \code{assembly = "hg19"}.
#' @param style Character value describing the style of arches.
#' Default value is \code{style = "2D"}. Options are:
#' \itemize{
#' \item{\code{"2D"}: }{Arches will be drawn in a 2-dimensional style.}
#' \item{\code{"3D"}: }{Arches will be drawn in a 3-dimensional style.}
#' }
#' @param flip Logical value indicating whether to reflect arches over
#' the x-axis. Default value is \code{flip = FALSE}.
#' @param curvature Numeric indicating the number of points along the
#' arch curvature. Default value is \code{curvature = 5}.
#' @param archHeight Single numeric value or numeric vector specifying
#' the arch heights. When NULL, all arches will be the same height,
#' filling up the given plot area
#' @param fill Character value(s) as a single value, vector, or palette
#' specifying fill colors of arches. Default value is \code{fill = #1f4297"}.
#' @param colorby A "\link[BentoBox]{colorby}" object specifying
#' information for scaling colors in data.
#' @param linecolor A character value specifying the color of the lines
#' outlining arches. Default value is \code{linecolor = NA}.
#' Special options include:
#' \itemize{
#' \item{\code{NA}: }{No line color.}
#' \item{\code{"fill"}: }{Same color as \code{fill}.}
#' }
#' @param alpha Numeric value specifying transparency.
#' Default value is \code{alpha = 0.4}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param clip A logical value indicating whether to clip any
#' arches that get cutoff in the given genomic region.
#' Default value is \code{clip = FALSE}.
#' @param baseline Logical value indicating whether to include
#' a baseline along the x-axis. Default value is \code{baseline = FALSE}.
#' @param baseline.color Baseline color.
#' Default value is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width.
#' Default value is \code{baseline.lwd = 1}.
#' @param x A numeric or unit object specifying BEDPE arches plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying BEDPE arches plot y-location.
#' The character value will
#' place the BEDPE arches plot y relative to the bottom of the most
#' recently plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying BEDPE arches plot width.
#' @param height A numeric or unit object specifying BEDPE arches plot height.
#' @param just Justification of BEDPE arches plot relative to its (x, y)
#' location. If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if
#' \code{x}, \code{y}, \code{width}, or \code{height} are only given as
#' numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be
#' produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_params} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_arches} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load BEDPE data
#' data("bb_bedpeData")
#'
#' ## Set the coordinates
#' params = bb_params(chrom = "chr21",
#'                    chromstart = 27900000, chromend = 30700000,
#'                    width = 7)
#'
#' ## Create a page
#' bb_pageCreate(width = 7.5, height = 2.1, default.units = "inches")
#'
#' ## Add a length column to color by
#' bb_bedpeData$length = (bb_bedpeData$start2 - bb_bedpeData$start1) / 1000
#'
#' ## Translate lengths into heights
#' heights = bb_bedpeData$length / max(bb_bedpeData$length)
#'
#' ## Plot the data
#' bedpePlot <- bb_plotBedpeArches(data = bb_bedpeData, params = params,
#'                                 fill = colorRampPalette(c("dodgerblue2", "firebrick2")),
#'                                 linecolor = "fill",
#'                                 colorby = colorby("length"),
#'                                 archHeight = heights, alpha = 1,
#'                                 x = 0.25, y = 0.25,  height = 1.5,
#'                                 just = c("left", "top"),
#'                                 default.units = "inches")
#'
#' ## Annotate genome label
#' bb_annoGenomeLabel(plot = bedpePlot, x = 0.25, y = 1.78, scale = "Mb")
#'
#'
#' ## Annotate heatmap legend
#' bb_annoHeatmapLegend(plot = bedpePlot, fontcolor = "black",
#'                      x = 7.0, y = 0.25,
#'                      width = 0.10, height = 1, fontsize = 10)
#'
#' ## Add the heatmap legend title
#' bb_plotText(label = "Kb", rot = 90, x = 6.9, y = 0.75,
#'             just = c("center","center"),
#'             fontsize = 10)
#'
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @details
#' A BEDPE Arches plot can be placed on a BentoBox coordinate page by providing plot placement parameters:
#' \preformatted{
#' bb_plotBedpeArches(data chrom,
#'                    chromstart = NULL, chromend = NULL,
#'                    x, y, width, height, just = c("left", "top"),
#'                    default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated BEDPE Arches plot by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotBedpeArches(data, chrom,
#'                    chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
bb_plotBedpeArches <- function(data, chrom, chromstart = NULL, chromend = NULL,
                               assembly = "hg19", style = "2D", flip = FALSE,
                               curvature = 5, archHeight = NULL,
                               fill = "#1f4297", colorby = NULL,
                               linecolor = NA, alpha = 0.4, bg = NA,
                               clip = FALSE, baseline = FALSE,
                               baseline.color = "grey", baseline.lwd = 1,
                               x = NULL, y = NULL, width = NULL, height = NULL,
                               just = c("left", "top"),
                               default.units = "inches", draw = TRUE,
                               params = NULL, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors
  errorcheck_bbArches <- function(bedpe, arches_plot, style, colorby){
    ## Can't have only one NULL chromstart or chromend
    if ((is.null(arches_plot$chromstart) & !is.null(arches_plot$chromend)) |
        (is.null(arches_plot$chromend) & !is.null(arches_plot$chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.", call. = FALSE)

    }


    if (!is.null(arches_plot$chromstart) & !is.null(arches_plot$chromend)){

      if (arches_plot$chromstart == arches_plot$chromend){
        stop("Genomic region is 0 bp long.", call. = FALSE)
      }

      ## chromend > chromstart
      if (arches_plot$chromend < arches_plot$chromstart){

        stop("\'chromstart\' should not be larger than \'chromend\'.",
             call. = FALSE)


      }

    }

    if (!style %in% c("3D", "2D")){
      stop("Invalid \'style\' input. Options are \'3D\' and \'2D\'.",
           call. = FALSE)
    }

    if (!is.null(colorby)){

      if (!any(colnames(bedpe) == colorby$column)){
        stop("Colorby column not found in data. Check colorby column name.",
             call. = FALSE)
      }

      if (length(which(colnames(bedpe) == colorby$column)) > 1){
        stop("Multiple matching colorby columns found in data. Please provide colorby column name with only one occurrence.", call. = FALSE)
      }
    }
}

  ## Define a function that will produce a yscale for arches based on height
  height_yscale <- function(heights, flip){

    if (length(heights) == 1){

      if (flip == FALSE){
        yscale <- c(0, heights)
      } else {
        yscale <- c(heights, 0)
      }

    } else {

      if (flip == FALSE){
        yscale <- c(0, max(heights))
      } else {
        yscale <- c(max(heights), 0)
      }


    }
    return(yscale)
  }

  ## Define a function that normalizes arch heights
  normHeights <- function(height, min, max){

    ## First normalize from 0 to 1
    newHeight <- (height - min)/(max-min)

    ## Then scale to a range of 1.38
    finalHeight <- newHeight * 1.38

    return(finalHeight)

  }

  ## Define a function that creates ribbon arch grobs
  drawRibbons <- function(df, style, arch, flip, transp, gp){

    x1 <- df[[2]]
    x2 <- df[[3]]
    y1 <- df[[5]]
    y2 <- df[[6]]
    fillCol <- df[[7]]
    lineCol <- df$linecolor
    outerHeight <- as.numeric(df[[length(df)]])
    innerHeight <- outerHeight - 0.01
    gp$fill <- fillCol
    gp$col <- lineCol
    gp$alpha <- transp

    if (style == "3D"){
      x1 <- df[[3]]
      x2 <- df[[2]]
    }

    ## Designate bezier control points
    # innerX = unit(seq(x2, y1, length.out = arch)[c(1:2, (arch-1):arch)],
    #               "native")
    innerX = unit(seq(x2, y1, length.out = arch)[c(1,2, seq((arch-1),arch))],
                  "native")
    # outerX = unit(seq(x1, y2, length.out = arch)[c(1:2, (arch-1):arch)],
    #               "native")
    outerX = unit(seq(x1, y2, length.out = arch)[c(1,2, seq((arch-1),arch))],
                  "native")

    if (flip == FALSE){
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
    archGrob <- polygonGrob(x = unit(c(convertX(outerBP$x, "native"),
                                       rev(convertX(innerBP$x, "native"))),
                                     "native"),
                            y = unit(c(convertY(outerBP$y, "npc"),
                                       rev(convertY(innerBP$y, "npc"))),
                                     "npc"),
                            gp = gp)
    assign("arches_grobs",
           addGrob(gTree = get("arches_grobs", envir = bbEnv),
                   child = archGrob), envir = bbEnv)


  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(clip)) clip <- NULL
  if(missing(style)) style <- NULL
  if(missing(curvature)) curvature <- NULL
  if(missing(flip)) flip <- NULL
  if(missing(fill)) fill <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(alpha)) alpha <- NULL
  if(missing(baseline)) baseline <- NULL
  if(missing(baseline.color)) baseline.color <- NULL
  if(missing(baseline.lwd)) baseline.lwd <- NULL
  if(missing(bg)) bg <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if bedpe/chrom arguments are missing (could be in object)
  if(!hasArg(data)) data <- NULL
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_archInternal <- structure(list(data = data, chrom = chrom,
                                    chromstart = chromstart,
                                    chromend = chromend,
                                    clip = clip, archHeight = archHeight,
                                    style = style,
                                    curvature = curvature, flip = flip,
                                    fill = fill, linecolor = linecolor,
                                    colorby = colorby, assembly = assembly,
                                    alpha = alpha, baseline = baseline,
                                    baseline.color = baseline.color,
                                    baseline.lwd = baseline.lwd,
                                    bg = bg, x = x, y = y, width = width,
                                    height = height, just = just,
                                    default.units = default.units,
                                    draw = draw, gp = gpar()),
                               class = "bb_archInternal")

  bb_archInternal <- parseParams(bb_params = params,
                                 object_params = bb_archInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_archInternal$clip)) bb_archInternal$clip <- FALSE
  if(is.null(bb_archInternal$style)) bb_archInternal$style <- "2D"
  if(is.null(bb_archInternal$curvature)) bb_archInternal$curvature <- 5
  if(is.null(bb_archInternal$flip)) bb_archInternal$flip <- FALSE
  if(is.null(bb_archInternal$fill)) bb_archInternal$fill <- "#1f4297"
  if(is.null(bb_archInternal$linecolor)) bb_archInternal$linecolor <- NA
  if(is.null(bb_archInternal$assembly)) bb_archInternal$assembly <- "hg19"
  if(is.null(bb_archInternal$alpha)) bb_archInternal$alpha <- 0.4
  if(is.null(bb_archInternal$baseline)) bb_archInternal$baseline <- FALSE
  if(is.null(bb_archInternal$baseline.color)) bb_archInternal$baseline.color <- "grey"
  if(is.null(bb_archInternal$baseline.lwd)) bb_archInternal$baseline.lwd <- 1
  if(is.null(bb_archInternal$bg)) bb_archInternal$bg <- NA
  if(is.null(bb_archInternal$just)) bb_archInternal$just <- c("left", "top")
  if(is.null(bb_archInternal$default.units)) bb_archInternal$default.units <- "inches"
  if(is.null(bb_archInternal$draw)) bb_archInternal$draw <- TRUE

  ## Set gp
  bb_archInternal$gp <- setGP(gpList = bb_archInternal$gp,
                              params = bb_archInternal, ...)

  # ======================================================================================================================================================================================
  # CHECK ARGUMENT ERRORS
  # ======================================================================================================================================================================================
  if (!is.null(bb_archInternal$colorby)){
    if(class(bb_archInternal$colorby) != "bb_colorby"){
      stop("\"colorby\" not of class \"bb_colorby\". Input colorby information with \"colorby()\".", call. = FALSE)
    }
  }

  if(is.null(bb_archInternal$data)) stop("argument \"data\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_archInternal$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  arches_plot <- structure(list(bedpe = NULL, chrom = bb_archInternal$chrom,
                                chromstart = bb_archInternal$chromstart,
                                chromend = bb_archInternal$chromend,
                                assembly = bb_archInternal$assembly,
                                color_palette = NULL,
                                zrange = bb_archInternal$colorby$range,
                                x = bb_archInternal$x, y = bb_archInternal$y,
                                width = bb_archInternal$width,
                                height = bb_archInternal$height,
                                just = bb_archInternal$just, grobs = NULL),
                           class = "bb_arches")
  attr(x = arches_plot, which = "plotted") <- bb_archInternal$draw

  # ======================================================================================================================================================================================
  # CHECK PLACEMENT ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = arches_plot)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  arches_plot$assembly <- parse_bbAssembly(assembly = arches_plot$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  arches_plot <- defaultUnits(object = arches_plot,
                              default.units = bb_archInternal$default.units)

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if ("data.frame" %in% class(bb_archInternal$data)){
    bedpe <- as.data.frame(bb_archInternal$data)
  } else {
    bedpe <- as.data.frame(data.table::fread(bb_archInternal$data))
  }

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bbArches(bedpe = bedpe, arches_plot = arches_plot,
                      style = bb_archInternal$style,
                      colorby = bb_archInternal$colorby)

  # ======================================================================================================================================================================================
  # WHOLE CHROM DATA AND XSCALE
  # ======================================================================================================================================================================================

  if (is.null(arches_plot$chromstart) & is.null(arches_plot$chromend)){

    if (class(arches_plot$assembly$TxDb) == "TxDb"){
      txdbChecks <- TRUE
    } else {
      txdbChecks <- check_loadedPackage(package = arches_plot$assembly$TxDb,
                                        message = paste(paste0("`",
                                                               arches_plot$assembly$TxDb,
                                                               "`"),
                                                        "not loaded. Please install and load to plot full chromosome ribbon arches."))
    }

    xscale <- c(0, 1)
    if (txdbChecks == TRUE){
      if (class(arches_plot$assembly$TxDb) == "TxDb"){
        tx_db <- arches_plot$assembly$TxDb
      } else {
        tx_db <- eval(parse(text = arches_plot$assembly$TxDb))
      }

      assembly_data <- GenomeInfoDb::seqlengths(tx_db)

      if (!arches_plot$chrom %in% names(assembly_data)){
        txdbChecks <- FALSE
        warning(paste("Chromosome", paste0("'",
                                           arches_plot$chrom,
                                           "'"),
                      "not found in", paste0("`",
                                             arches_plot$assembly$TxDb,
                                             "`"),
                      "and data for entire chromosome cannot be plotted."),
                call. = FALSE)
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

    if (bb_archInternal$clip == TRUE){
      bedpe <- bedpe[which(bedpe[,1] == arches_plot$chrom
                           & bedpe[,4] == arches_plot$chrom
                           & bedpe[,2] >= arches_plot$chromstart
                           & bedpe[,3] <= arches_plot$chromend
                           & bedpe[,5] >= arches_plot$chromstart
                           & bedpe[,6] <= arches_plot$chromend),]
    } else {

      bedpe <- bedpe[which(bedpe[,1] == arches_plot$chrom
                           & bedpe[,4] == arches_plot$chrom
                           & ((bedpe[,3] >= arches_plot$chromstart
                               & bedpe[,3] <= arches_plot$chromend)
                              | (bedpe[,5] <= arches_plot$chromstart
                                 & bedpe[,5] >= arches_plot$chromend))),]
    }

  } else {

    bedpe <- data.frame(matrix(nrow = 0, ncol = 6))
  }

  arches_plot$bedpe <- bedpe

  # ======================================================================================================================================================================================
  # COLORBY AND COLORS
  # ======================================================================================================================================================================================

   if (!is.null(bb_archInternal$colorby) & nrow(bedpe) > 0){
    colorbyCol <- which(colnames(bedpe) == bb_archInternal$colorby$column)
    colorbyCol <- bedpe[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" & class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      colorbyCol <- as.numeric(colorbyCol)
    }

    if (is.null(bb_archInternal$colorby$range)){
      colorbyrange <- c(min(colorbyCol), max(colorbyCol))
      arches_plot$zrange <- colorbyrange
    }

    if (class(bb_archInternal$fill) == "function"){
      colors <- bb_maptocolors(colorbyCol, bb_archInternal$fill,
                               range = arches_plot$zrange)
      arches_plot$color_palette <- bb_archInternal$fill

    } else {
      colorbyColfac <- factor(colorbyCol)
      mappedColors <- rep(bb_archInternal$fill,
                          ceiling(length(levels(colorbyColfac))/length(bb_archInternal$fill)))

      names(mappedColors) <- levels(colorbyColfac)
      colors <- mappedColors[colorbyColfac]

    }

  } else {
    if (class(bb_archInternal$fill) == "function"){
      colors <- bb_archInternal$fill(nrow(bedpe))
    } else {
      if (length(bb_archInternal$fill) == 1){
        colors <- rep(bb_archInternal$fill, nrow(bedpe))
      } else {
        colors <- rep(bb_archInternal$fill,
                      ceiling(nrow(bedpe)/length(bb_archInternal$fill)))[seq(1, nrow(bedpe))]
      }
    }

  }

  bedpe <- bedpe[,c(seq(1, 6))]
  bedpe$color <- colors

  # Set actual line color to fill color if requested by user
  actuallinecolor = bb_archInternal$linecolor
  if (is.na(bb_archInternal$linecolor) == FALSE)
  {
    if (bb_archInternal$linecolor == "fill")
    {
      actuallinecolor = bedpe$color
    }
  }
  if (is.null(bb_archInternal$linecolor) == TRUE)
  {
    actuallinecolor = NA
  }
  bedpe$linecolor = actuallinecolor

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_arches",
                    length(grep(pattern = "bb_arches",
                                x = currentViewports)) + 1)

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

    add_bbViewport(vp_name)

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

  if (nrow(bedpe) > 0){
  if (is.null(bb_archInternal$archHeight)){
    bedpe$height <- rep(1, nrow(bedpe))
  } else if (length(bb_archInternal$archHeight) == 1){
    bedpe$height <- rep(bb_archInternal$archHeight, nrow(bedpe))
    yscale <- height_yscale(heights = bb_archInternal$archHeight,
                            flip = bb_archInternal$flip)
    vp$yscale <- yscale
  } else {
    bedpe$height <- bb_archInternal$archHeight[seq(1, nrow(bedpe))]
    yscale <- height_yscale(heights = bb_archInternal$archHeight,
                            flip = bb_archInternal$flip)
    vp$yscale <- yscale
  }

  if (length(bedpe$height) > 0){
    bedpe$normHeight <- lapply(bedpe$height, normHeights, min = 0,
                               max = max(bedpe$height))
  }


    }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
  # ======================================================================================================================================================================================

  backgroundGrob <- rectGrob(gp = gpar(fill = bb_archInternal$bg, col = NA),
                             name = "background")
  assign("arches_grobs", gTree(vp = vp, children = gList(backgroundGrob)),
         envir = bbEnv)

  # ======================================================================================================================================================================================
  # GROBS
  # ======================================================================================================================================================================================

  if (nrow(bedpe) > 0){

    if (bb_archInternal$baseline == TRUE){
      baselineGrob <- segmentsGrob(x0 = unit(0, "npc"), y0 = 0,
                                   x1 = unit(1, "npc"), y1 = 0,
                                   default.units = "native",
                                   gp = gpar(col = bb_archInternal$baseline.color,
                                             lwd = bb_archInternal$baseline.lwd))
      assign("arches_grobs",
             addGrob(gTree = get("arches_grobs", envir = bbEnv),
                     child = baselineGrob), envir = bbEnv)

    }

    invisible(apply(bedpe, 1, drawRibbons, style = bb_archInternal$style,
                    arch = bb_archInternal$curvature,
                    flip = bb_archInternal$flip,
                    transp = bb_archInternal$alpha, gp = bb_archInternal$gp))

  } else {

    if (txdbChecks == TRUE){
      warning("BEDPE data contains no values.", call. = FALSE)
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

  message(paste0("bb_arches[", vp$name, "]"))
  invisible(arches_plot)

}




