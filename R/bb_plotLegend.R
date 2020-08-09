#' plots a horizontal legend
#'
#' @param legend single value or vector of legend text
#' @param fillcolor if specified, this argument will produce boxes filled with the specified colors to appear beside the legend text
#' @param pch if specified, this argument will produce symbols to appear beside the legend text
#' @param lty if specified, this argument will produce line types to appear beside the legend text
#' @param title legend title
#' @param border option to add border around legend
#' @param x A numeric vector or unit object specifying x-location
#' @param y A numeric vector or unit object specifying y-location
#' @param width A numeric vector or unit object specifying width
#' @param height A numeric vector or unit object specifying height
#' @param just justification of scale viewport
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#' @param draw A logical value indicating whether graphics output should be produced

#' @export

bb_plotLegend <- function(legend, fillcolor = NULL, pch = NULL, lty = NULL, title = NULL, fontsize = 10, border = TRUE, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "bottom"),
                          default.units = "inches", draw = TRUE, ...){

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  legend_plot <- structure(list(width = width, height = height, x = x, y = y, justification = just, grobs = NULL), class = "bb_legend")
  attr(x = legend_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = legend_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  legend_plot <- defaultUnits(object = legend_plot, default.units = default.units)

  textHeight <- heightDetails(textGrob(label = "A", gp = gpar(fontsize = fontsize)))
  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_legend", length(grep(pattern = "bb_legend", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(0.125, "snpc"), width = unit(0.20, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   just = "center",
                   yscale = c(0, 0.125),
                   xscale = c(0, 0.20),
                   name = vp_name)
    pushViewport(vp)
    textHeight <- convertHeight(textHeight, unitTo = "native", valueOnly = T)
    upViewport()

    height <- 0.125

    if (draw == TRUE){

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
                   just = just,
                   yscale = c(0, height),
                   xscale = c(0, width),
                   name = vp_name)

    textHeight <- convertHeight(textHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("legend_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # GROBS
  # ======================================================================================================================================================================================

  ## Title
  if (!is.null(title)){
    remainingSpace <- height - textHeight*(length(legend)+1)
    spaceNo <- length(legend) + 2

    titleGrob <- textGrob(label = title, x = unit(0.5, "npc"), y = height - remainingSpace/spaceNo,
                          just = "top", gp = gpar(fontsize = fontsize), default.units = "native")
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = titleGrob), envir = bbEnv)

  } else {
    remainingSpace <- height - textHeight*length(legend)
    spaceNo <- length(legend) + 1
  }

  spaceHeight <- remainingSpace/spaceNo

  ycoords <- seq(from = spaceHeight, to = height, by = (spaceHeight + textHeight))[1:length(legend)]
  xcoord <- textHeight*2

  ## Colors
  ## Only take the first number of colors for the length of legend
  fillcolors <- fillcolor[1:length(legend)]

  ## Symbols
  if (!is.null(pch)){
    ## Only take the first number of symbols for the length of legend
    pchs <- pch[1:length(legend)]

    labelSymbols <- pointsGrob(x = textHeight, y = ycoords,
                               pch = pchs, size = unit(fontsize, "points"),
                               gp = gpar(fill = fillcolors), default.units = "native")
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = labelSymbols), envir = bbEnv)

    ## Lines
  } else if (!is.null(lty)) {
    ## Only take the first number of symbols for the length of legend
    ltys <- lty[1:length(legend)]
    xcoord <- 3*textHeight
    labelLines <- segmentsGrob(x0 = 0.5*textHeight, y0 = ycoords + 0.5*textHeight,
                               x1 = 2.5*textHeight, y1 = ycoords + 0.5*textHeight,
                               gp = gpar(lty = ltys, col = fillcolors), default.units = "native")
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = labelLines), envir = bbEnv)

  } else {

    labelColors <- rectGrob(x = 0.5*textHeight, y = ycoords,
                            width = textHeight, height = textHeight,
                            just = c("left", "bottom"), gp = gpar(fill = fillcolors, col = NA), default.units = "native")
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = labelColors), envir = bbEnv)
  }


  ## Text labels
  labelText <- textGrob(label = rev(legend), x = xcoord, y = ycoords + 0.5*textHeight,
                        just = "left", gp = gpar(fontsize = fontsize, ...), default.units = "native")
  assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = labelText), envir = bbEnv)


  ## Border
  if (border == TRUE){
    border <- rectGrob(gp = gpar(fill = NA, col = "black"))
    assign("legend_grobs", addGrob(get("legend_grobs", envir = bbEnv), child = border), envir = bbEnv)
  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

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
