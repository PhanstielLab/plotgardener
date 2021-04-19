#' Annotates a highlight box around a specified genomic region of a BentoBox plot
#'
#' @param plot Input BentoBox plot on which to annotate genomic region.
#' @param chrom Chromosome of region to be highlighted, as a string.
#' @param chromstart Integer start position on chromosome to be highlighted.
#' @param chromend Integer end position on chromosome to be highlighted.
#' @param fill A character value specifying highlight box fill color. Default value is \code{fill = "grey"}.
#' @param linecolor A character value specifying highlight box line color. Default value is \code{linecolor = NA}.
#' @param alpha Numeric value specifying color transparency. Default value is \code{alpha = 0.4}.
#' @param y A numeric, unit object, or character containing a "b" combined with a numeric value specifying square highlight box y-location. The character value will
#' place the highlight box y relative to the bottom of the most recently plotted BentoBox plot according to the units of the BentoBox page.
#' @param height A numeric or unit object specifying highlight box height.
#' @param just Justification of highlight box relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{y} or \code{height} are only given as numerics or numeric vectors. Default value is \code{default.units = "inches"}.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_highlight} object containing relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' bb_pageCreate(width = 7.5, height = 1.5, default.units = "inches")
#'
#' ## Plot and place a signal plot
#' data("bb_imrH3K27acData")
#' region <- bb_params(chrom = "chr21", chromstart = 28000000, chromend = 30300000,
#'                     range = c(0, 45))
#' signalPlot <- bb_plotSignal(data = bb_imrH3K27acData, params = region,
#'                             x = 0.5, y = 0.25, width = 6.5, height = 0.65,
#'                             just = c("left", "top"), default.units = "inches")
#'
#' ## Highlight genomic region on signal plot
#' bb_annoHighlight(plot = signalPlot,
#'                  chrom = "chr21", chromstart = 29000000, chromend = 29125000,
#'                  y = 0.25, height = 1, just = c("left", "top"), default.units = "inches")
#'
#' ## Plot text label
#' bb_plotText(label = "region of interest", fontsize = 8, fontcolor = "black",
#'             x = 3.5, y = 0.2, just = "bottom", default.units = "inches")
#'
#' ## Plot genome label
#' bb_plotGenomeLabel(chrom = "chr21", chromstart = 28000000, chromend = 30300000,
#'                    x = 0.5, y = 1.3, length = 6.5, default.units = "inches")
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @export
bb_annoHighlight <- function(plot, chrom, chromstart = NULL, chromend = NULL, fill = "grey", linecolor = NA, alpha = 0.4,
                             y, height, just = c("left", "top"), default.units = "inches", params = NULL, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_annoHighlight
  bb_errorcheck_annoHighlight <- function(object){

    if (!object$chrom %in% object$plot$chrom){
      stop(paste(object$chrom, "not found in input plot. Cannot annotate highlight."), call. = FALSE)
    }

    if(length(object$chrom) > 1){
      stop("Cannot highlight multiple chromosome regions in one function call.", call. = FALSE)
    }

  }

  ## Define a function that converts chromstart/chromend to page units for multi-chromosome manhattan plot
  convertManhattan <- function(object, manhattanPlot){

    ## Get assembly data
    if (class(object$assembly$TxDb) == "TxDb"){
      txdbChecks <- TRUE
    } else {
      txdbChecks <- suppressWarnings(check_loadedPackage(package = object$assembly$TxDb, message = NULL))
    }


    if (txdbChecks == TRUE){
      if (class(object$assembly$TxDb) == "TxDb"){
        tx_db <- object$assembly$TxDb
      } else {
        tx_db <- eval(parse(text = object$assembly$TxDb))
      }

      assembly_data <- as.data.frame(setDT(as.data.frame(seqlengths(tx_db)), keep.rownames = TRUE))
      assembly_data <- assembly_data[which(assembly_data[,1] %in% manhattanPlot$chrom),]

      ## Get the offset based on spacer for the assembly
      offsetAssembly <- spaceChroms(assemblyData = assembly_data, space = manhattanPlot$space)
      offsetAssembly <- offsetAssembly[which(offsetAssembly$chrom == object$chrom),]

      ## Convert chromstart and chromend to chrom offsetAssembly range
      oldRange <- offsetAssembly[,2] - 1
      newRange <- offsetAssembly[,4] - offsetAssembly[,3]
      newStart <- (((object$chromstart - 1) * newRange) / oldRange) + offsetAssembly[,3]
      newEnd <- (((object$chromend - 1) * newRange) / oldRange) + offsetAssembly[,3]

      ## Convert new chromstart and chromend to page units
      start <- convertX(unit(newStart, "native"), unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)
      end <- convertX(unit(newEnd, "native"), unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)

      return(list(start, end))

    }

  }

  ## Define a function that resets y to a top-based justification
  y_Pagetop <- function(y, height, just){
    page_height <- get("page_height", envir = bbEnv)
    page_units <- get("page_units", envir = bbEnv)
    y <- convertY(unit(page_height, units = page_units) - convertY(y, unitTo = page_units), unitTo = page_units)
    height <- convertHeight(height, unitTo = page_units)

    if (length(just == 2)){

      if ("left" %in% just & "center" %in% just){
        topy <- y + 0.5*height
      } else if ("right" %in% just & "center" %in% just){
        topy <- y + 0.5*height
      } else if ("center" %in% just & "bottom" %in% just){
        topy <- y + height
      } else if ("center" %in% just & "top" %in% just){
        topy <- y
      } else if ("left" %in% just & "top" %in% just){
        topy <- y
      } else if ("right" %in% just & "top" %in% just){
        topy <- y
      } else if ("left" %in% just & "bottom" %in% just){
        topy <- y + height
      } else if ("right" %in% just & "bottom" %in% just){
        topy <- y + height
      } else {
        topy <- y + 0.5*height
      }

    } else if (length(just == 1)){

      if (just == "left"){
        topy <- y + 0.5*height
      } else if (just == "right"){
        topy <- y + 0.5*height
      } else if (just == "bottom"){
        topy <- y + height
      } else if (just == "top"){
        topy <- y
      } else {
        topy <- y + 0.5*height
      }

    }

    return(topy)

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(fill)) fill <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(alpha)) alpha <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if arguments are missing (could be in object)
  if(!hasArg(plot)) plot <- NULL
  if(!hasArg(chrom)) chrom <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(height)) height <- NULL

  ## Compile all parameters into an internal object
  bb_highlightInternal <- structure(list(plot = plot, chrom = chrom, chromstart = chromstart, chromend = chromend,
                                         fill = fill, linecolor = linecolor, alpha = alpha, y = y, height = height,
                                         just = just, default.units = default.units), class = "bb_highlightInternal")

  bb_highlightInternal <- parseParams(bb_params = params, object_params = bb_highlightInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_highlightInternal$fill)) bb_highlightInternal$fill <- "grey"
  if(is.null(bb_highlightInternal$linecolor)) bb_highlightInternal$linecolor <- NA
  if(is.null(bb_highlightInternal$alpha)) bb_highlightInternal$alpha <- 0.4
  if(is.null(bb_highlightInternal$just)) bb_highlightInternal$just <- c("left", "top")
  if(is.null(bb_highlightInternal$default.units)) bb_highlightInternal$default.units <- "inches"

  ## Set gp
  bb_highlightInternal$gp <- gpar(fill = bb_highlightInternal$fill, col = bb_highlightInternal$linecolor,
                                  alpha = bb_highlightInternal$alpha)
  bb_highlightInternal$gp <- setGP(gpList = bb_highlightInternal$gp, params = bb_highlightInternal, ...)

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_highlight <- structure(list(chrom = bb_highlightInternal$chrom, chromstart = bb_highlightInternal$chromstart, chromend = bb_highlightInternal$chromend,
                                 assembly = bb_highlightInternal$plot$assembly, x = NULL, y = bb_highlightInternal$y, width = NULL, height = bb_highlightInternal$height,
                                 just = bb_highlightInternal$just, grobs = NULL), class = "bb_highlight")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot add highlight annotation without a BentoBox page.")
  if(is.null(bb_highlightInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_highlightInternal$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_highlight$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_highlight$height)) stop("argument \"height\" is missing, with no default.", call. = FALSE)

  bb_errorcheck_annoHighlight(object = bb_highlightInternal)
  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  page_units <- get("page_units", envir = bbEnv)

  if (!"unit" %in% class(bb_highlight$y)){

    ## Check for "below" y-coords
    if (all(grepl("b", bb_highlight$y)) == TRUE){
      if (any(grepl("^[ac-zA-Z]+$", bb_highlight$y)) == TRUE){
        stop("\'below\' y-coordinate(s) detected with additional letters. Cannot parse y-coordinate(s).", call. = FALSE)
      }

      if(any(is.na(as.numeric(gsub("b","", bb_highlight$y))))){
        stop("\'below\' y-coordinate(s) does not have a numeric associated with it. Cannot parse y-coordinate(s).", call. = FALSE)
      }

      bb_highlight$y <- unit(unlist(lapply(bb_highlight$y, plot_belowY)), page_units)

    } else {

      if (!is.numeric(bb_highlight$y)){

        stop("y-coordinate is neither a unit object or a numeric value. Cannot plot text.", call. = FALSE)

      }

      if (is.null(bb_highlightInternal$default.units)){

        stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      bb_highlight$y <- unit(bb_highlight$y, bb_highlightInternal$default.units)

    }

  }

  if (!"unit" %in% class(bb_highlight$height)){

    if (!is.numeric(bb_highlight$height)){

      stop("Height is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_highlightInternal$default.units)){

      stop("Height detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_highlight$height <- unit(bb_highlight$height, bb_highlightInternal$default.units)

  }

  if (class(bb_highlightInternal$plot) == "bb_genes"){

    plotVP <- bb_highlightInternal$plot$grobs$children$background$vp

  } else if (class(bb_highlightInternal$plot) == "bb_hicTriangle" | class(bb_highlightInternal$plot) == "bb_hicRectangle"){

    plotVP <- bb_highlightInternal$plot$outsideVP

  } else {

    plotVP <- bb_highlightInternal$plot$grobs$vp

  }

  # ======================================================================================================================================================================================
  # WHOLE CHROM
  # ======================================================================================================================================================================================

  if (is.null(bb_highlight$chromstart) & is.null(bb_highlight$chromend)){

    if (class(bb_highlight$assembly$TxDb) == "TxDb"){
      txdbChecks <- TRUE
    } else {
      txdbChecks <- check_loadedPackage(package = bb_highlight$assembly$TxDb, message = paste(paste0("`", bb_highlight$assembly$TxDb,"`"),
                                                                                              "not loaded. Please install and load to annotate full chromosome region of plot."))
    }

    if (txdbChecks == TRUE){

      if (class(bb_highlight$assembly$TxDb) == "TxDb"){
        tx_db <- bb_highlight$assembly$TxDb
      } else {
        tx_db <- eval(parse(text = bb_highlight$assembly$TxDb))
      }

      assembly_data <- seqlengths(tx_db)
      if (!bb_highlight$chrom %in% names(assembly_data)){
        warning(paste("Chromosome", paste0("'", bb_highlight$chrom, "'"), "not found in", paste0("`", bb_highlight$assembly$TxDb, "`"), "and data for entire chromosome cannot be highlighted."), call. = FALSE)
      } else {
        bb_highlight$chromstart <- 1
        bb_highlight$chromend <- assembly_data[[bb_highlight$chrom]]
      }

    }

  }

  # ======================================================================================================================================================================================
  # DETERMINE X AND WIDTH BASED ON GENOMIC REGION
  # ======================================================================================================================================================================================

  ## Get plot viewport
  if (class(bb_highlightInternal$plot) == "bb_genes"){

    plotVP <- bb_highlightInternal$plot$grobs$children$background$vp

  } else if (class(bb_highlightInternal$plot) == "bb_hicTriangle" | class(bb_highlightInternal$plot) == "bb_hicRectangle"){

    plotVP <- bb_highlightInternal$plot$outsideVP

  } else {

    plotVP <- bb_highlightInternal$plot$grobs$vp

  }

  ## Convert plot viewport to bottom left to get left position on entire page
  plotVP_bottomLeft <- vp_bottomLeft(viewport = plotVP)

  ## Convert plot genomic coordinates to position on page
  seekViewport(plotVP$name)

  if (class(bb_highlightInternal$plot) == "bb_manhattan"){

    ## Multiple chromosome manhattan plot
    if (length(bb_highlightInternal$plot$chrom) > 1){

      convertedCoords <- convertManhattan(object = bb_highlight, manhattanPlot = bb_highlightInternal$plot)
      start <- convertedCoords[[1]]
      end <- convertedCoords[[2]]

    } else {
      start <- convertX(unit(bb_highlight$chromstart, "native"), unitTo = page_units, valueOnly = TRUE)
      end <- convertX(unit(bb_highlight$chromend, "native"), unitTo = page_units, valueOnly = TRUE)
    }

  } else {
    start <- convertX(unit(bb_highlight$chromstart, "native"), unitTo = page_units, valueOnly = TRUE)
    end <- convertX(unit(bb_highlight$chromend, "native"), unitTo = page_units, valueOnly = TRUE)
  }

  width <- end - start
  seekViewport("bb_page")

  ## Add additional page units to start
  start <- as.numeric(plotVP_bottomLeft[[1]]) + start

  # ======================================================================================================================================================================================
  # USE JUSTIFICATION TO DETERMINE PAGE-ADJUSTED TOP Y-COORD
  # ======================================================================================================================================================================================

  top_y <- y_Pagetop(y = bb_highlight$y, height = bb_highlight$height, just = bb_highlight$just)

  # ======================================================================================================================================================================================
  # HIGHLIGHT GROB
  # ======================================================================================================================================================================================

  name <- paste0("bb_highlight", length(grep(pattern = "bb_highlight", x = grid.ls(print = FALSE, recursive = FALSE))) + 1)
  highlightGrob <- grid.rect(x = unit(start, page_units), y = top_y, width = unit(width, page_units), height = bb_highlight$height,
                             just = c("left", "top"),
                             gp = bb_highlightInternal$gp, name = name)

  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================

  bb_highlight$grobs <- highlightGrob

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_highlight[", name, "]"))
  invisible(bb_highlight)

}
