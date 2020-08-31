#' plots a legend
#'
#' @param legend single value or vector of legend text#'
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param orientation "v" (vertical) or "h" (horizontal) orientation
#' @param fillcolor if specified, this argument will produce boxes filled with the specified colors to appear beside the legend text
#' @param pch if specified, this argument will produce symbols to appear beside the legend text
#' @param lty if specified, this argument will produce line types to appear beside the legend text
#' @param title legend title
#' @param border option to add border around legend
#' @param bg background color of legend
#' @param x A numeric vector or unit object specifying x-location
#' @param y A numeric vector or unit object specifying y-location
#' @param width A numeric vector or unit object specifying width
#' @param height A numeric vector or unit object specifying height
#' @param just justification of scale viewport
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#' @param draw A logical value indicating whether graphics output should be produced

#' @export

bb_plotLegend <- function(legend, params = NULL, orientation = "v", fillcolor = NULL, pch = NULL, lty = NULL, title = NULL, fontsize = 10, border = TRUE,
                          bg = NA, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, ...){

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(orientation)) orientation <- NULL
  if(missing(fontsize)) fontsize <- NULL
  if(missing(border)) border <- NULL
  if(missing(bg)) bg <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if legend argument is missing (could be in object)
  if(!hasArg(legend)) legend <- NULL

  ## Compile all parameters into an internal object
  bb_legInternal <- structure(list(legend = legend, orientation = orientation, fillcolor = fillcolor, pch = pch, lty = lty, title = title, fontsize = fontsize,
                                   border = border, bg = bg, x = x, y = y, width = width, height = height, just = just, default.units = default.units, draw = draw), class = "bb_legInternal")

  bb_legInternal <- parseParams(bb_params = params, object_params = bb_legInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_legInternal$orientation)) bb_legInternal$orientation <- "v"
  if(is.null(bb_legInternal$fontsize)) bb_legInternal$fontsize <- 10
  if(is.null(bb_legInternal$border)) bb_legInternal$border <- TRUE
  if(is.null(bb_legInternal$bg)) bb_legInternal$bg <- NA
  if(is.null(bb_legInternal$just)) bb_legInternal$just <- c("left", "top")
  if(is.null(bb_legInternal$default.units)) bb_legInternal$default.units <- "inches"
  if(is.null(bb_legInternal$draw)) bb_legInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  legend_plot <- structure(list(width = bb_legInternal$width, height = bb_legInternal$height, x = bb_legInternal$x, y = bb_legInternal$y,
                                justification = bb_legInternal$just, grobs = NULL), class = "bb_legend")
  attr(x = legend_plot, which = "plotted") <- bb_legInternal$draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================
  if (is.null(bb_legInternal$legend)) stop("argument \"legend\" is missing, with no default.", call. = FALSE)

  check_placement(object = legend_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  legend_plot <- defaultUnits(object = legend_plot, default.units = bb_legInternal$default.units)

  textHeight <- heightDetails(textGrob(label = "A", gp = gpar(fontsize = bb_legInternal$fontsize)))
  textGrobs <- lapply(bb_legInternal$legend, textGrob, gp = gpar(fontsize = bb_legInternal$fontsize))
  textWidths <- lapply(textGrobs, widthDetails)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_legend", length(grep(pattern = "bb_legend", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(legend_plot$x) & is.null(legend_plot$y)){

    vp <- viewport(height = unit(0.125, "snpc"), width = unit(0.20, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   just = "center",
                   yscale = c(0, 0.125),
                   xscale = c(0, 0.20),
                   name = vp_name)
    pushViewport(vp)
    textHeight <- convertHeight(textHeight, unitTo = "native", valueOnly = T)
    textWidths <- lapply(textWidths, convertWidth, unitTo = "native", valueOnly = T)
    upViewport()

    height <- 0.125
    width <- 0.20

    if (bb_legInternal$draw == TRUE){

      vp$name <- "bb_legend1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = legend_plot)
    height <- convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    width <- convertWidth(page_coords$width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   just = bb_legInternal$just,
                   yscale = c(0, height),
                   xscale = c(0, width),
                   name = vp_name)

    textHeight <- convertHeight(textHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    textWidths <- lapply(textWidths, convertWidth, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("legend_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # GROBS
  # ======================================================================================================================================================================================

  ## Border
  if (bb_legInternal$border == TRUE){
    border <- rectGrob(gp = gpar(fill = bb_legInternal$bg, col = "black"))
  } else {
    border <- rectGrob(gp = gpar(fill = bb_legInternal$bg, col = NA))
  }

  assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = border), envir = bbEnv)

  ## Title
  if (!is.null(bb_legInternal$title)){

    if (bb_legInternal$orientation == "h"){

      remainingSpace <- height - textHeight*2
      spaceNo <- 3

    } else {

      remainingSpace <- height - textHeight*(length(bb_legInternal$legend)+1)
      spaceNo <- length(bb_legInternal$legend) + 2

    }

    titleGrob <- textGrob(label = bb_legInternal$title, x = unit(0.5, "npc"), y = height - remainingSpace/spaceNo,
                          just = "top", gp = gpar(fontsize = bb_legInternal$fontsize), default.units = "native")
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = titleGrob), envir = bbEnv)

  } else {

    if (bb_legInternal$orientation == "h"){

      remainingSpace <- height - textHeight
      spaceNo <- 2

    } else {
      remainingSpace <- height - textHeight*length(bb_legInternal$legend)
      spaceNo <- length(legend) + 1

    }

  }


  ## Colors
  ## Only take the first number of colors for the length of legend
  fillcolors <- bb_legInternal$fillcolor[1:length(bb_legInternal$legend)]
  ## Spacing and label coordinates
  spaceHeight <- remainingSpace/spaceNo

  if (bb_legInternal$orientation == "h"){
    ycoords <- spaceHeight
    remainingWidth <- width - (sum(unlist(textWidths)) + textHeight*length(bb_legInternal$legend))
    width_spaces <- 2*length(bb_legInternal$legend) + 1
    widthSpace <- remainingWidth/width_spaces
    if (widthSpace < 0){
      widthSpace <- 0.1*textHeight
    }

    addTexts <- function(factor, textWidths){
      totalText <- sum(textWidths[0:factor])
      return(totalText)
    }

    textWidth <- c(0, unlist(textWidths)[1:(length(bb_legInternal$legend)-1)])
    df <- as.data.frame(cbind("factor" = 0:(length(bb_legInternal$legend)-1),
                "firstSpace" = rep(widthSpace, length(bb_legInternal$legend))))
    df$totalText <- unlist(lapply(df$factor, addTexts, textWidths = unlist(textWidths)))

    xcoords <- df$firstSpace + df$factor*(2*widthSpace + textHeight) + df$totalText

    xpch <- xcoords + 0.5*textHeight
    ypch <- rep(spaceHeight+0.5*textHeight, length(bb_legInternal$legend))

    x0s <- xcoords
    y0s <- spaceHeight + 0.5*textHeight
    x1s <- x0s + textHeight
    y1s <- spaceHeight + 0.5*textHeight


  } else {

    ycoords <- seq(from = spaceHeight, to = height, by = (spaceHeight + textHeight))[1:length(bb_legInternal$legend)]
    xcoords <- textHeight*2

    xpch <- rep(textHeight*2, length(bb_legInternal$legend))
    ypch <- ycoords + 0.5*textHeight

    x0s <- textHeight
    y0s <- ycoords + 0.5*textHeight
    x1s <- 3*textHeight
    y1s <- ycoords + 0.5*textHeight

    fillcolors <- rev(fillcolors)
  }




  ## Symbols
  if (!is.null(bb_legInternal$pch)){
    ## Only take the first number of symbols for the length of legend
    pchs <- bb_legInternal$pch[1:length(bb_legInternal$legend)]

    labelSymbols <- pointsGrob(x = xpch, y = ypch,
                               pch = pchs, size = unit(bb_legInternal$fontsize, "points"),
                               gp = gpar(fill = fillcolors), default.units = "native")
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = labelSymbols), envir = bbEnv)

    ## Lines
  } else if (!is.null(lty)) {
    ## Only take the first number of symbols for the length of legend
    ltys <- lty[1:length(legend)]
    labelLines <- segmentsGrob(x0 = x0s, y0 = y0s,
                               x1 = x1s, y1 = y1s,
                               gp = gpar(lty = ltys, col = fillcolors), default.units = "native")
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = labelLines), envir = bbEnv)

  } else {

    labelColors <- rectGrob(x = xcoords, y = ycoords,
                            width = textHeight, height = textHeight,
                            just = c("left", "bottom"), gp = gpar(fill = fillcolors, col = NA), default.units = "native")
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = labelColors), envir = bbEnv)
  }


  ## Text labels
  if (bb_legInternal$orientation == "h"){

    labelText <- textGrob(label = bb_legInternal$legend, x = xcoords + textHeight + widthSpace, y = ycoords + 0.5*textHeight,
                          just = "left", gp = gpar(fontsize = bb_legInternal$fontsize,...), default.units = "native")

  } else {

    labelText <- textGrob(label = rev(bb_legInternal$legend), x = xcoords + 2*textHeight, y = ycoords + 0.5*textHeight,
                          just = "left", gp = gpar(fontsize = bb_legInternal$fontsize,...), default.units = "native")
  }

  assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = labelText), envir = bbEnv)


  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_legInternal$draw == TRUE){

    grid.draw(get("legend_grobs", envir = bbEnv))

  }
  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  legend_plot$grobs <-  get("legend_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(legend_plot)

}
