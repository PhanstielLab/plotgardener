#' Add a color scale legend for heatmap-style plots
#'
#' @param plot Heatmap-style plot object to add heatmap legend for.
#' @param orientation A string specifying legend orientation. Default value is \code{orientation = "v"}. Options are:
#' \itemize{
#' \item{\code{"v"}: }{Vertical legend orientation.}
#' \item{\code{"h"}: }{Horizontal legend orientation.}
#' }
#' @param fontsize A numeric specifying text fontsize in points. Default value is \code{fontsize = 8}.
#' @param fontcolor Character value specfying text fontcolor. Default value is \code{fontcolor = "dark grey"}.
#' @param border Logical value indicating whether to add a border around heatmap legend. Default value is \code{border = FALSE}.
#' @param x A numeric or unit object specifying x-location of legend.
#' @param y A numeric or unit object specifying y-location of legend.
#' @param width A numeric or unit object specifying width of legend.
#' @param height A numeric or unit object specifying height of legend.
#' @param just Justification of heatmap legend relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification. Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{x}, \code{y}, \code{width}, or \code{height} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_heatmapLegend} object with relevant color value, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data
#' data("bb_hicData")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 4, height = 3.5, default.units = "inches", xgrid = 0, ygrid = 0)
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- bb_plotHicSquare(hicData = bb_hicData, resolution = 10000, zrange = c(0, 70),
#'                             chrom = "chr21", chromstart = 28000000, chromend = 30300000,
#'                             x = 0.5, y = 0.5, width = 2.5, height = 2.5,
#'                             just = c("left", "top"), default.units = "inches")
#'
#' ## Add heatmap legend
#' bb_annoHeatmapLegend(plot = hicPlot,
#'                      x = 3.2, y = 0.5, width = 0.25, height = 0.75,
#'                      just = c("left", "top"), default.units = "inches")
#'
#' @export
bb_annoHeatmapLegend <- function(plot, orientation = "v", fontsize = 8, fontcolor = "dark grey", border = FALSE, x, y, width, height,
                                 just = c("left", "top"), default.units = "inches", params = NULL, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_annoHeatmapLegend
  errorcheck_bb_annoHeatmapLegend <- function(bb_heatmapLegend, orientation){

    ## checking min_val and max val
    if (as.numeric(bb_heatmapLegend$min_val) > as.numeric(bb_heatmapLegend$max_val)){

      warning("/'min_val/' is larger than /'max_val/'. Legend labels may be incorrect.", call. = FALSE)

    }

    ## check proper orientation
    if (!orientation %in% c("v", "h")){

      stop("Invalid /'orientation/' parameter. Options are 'v' or 'h'.", call. = FALSE)

    }

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(border)) border <- NULL
  if(missing(orientation)) orientation <- NULL
  if(missing(fontsize)) fontsize <- NULL
  if(missing(fontcolor)) fontcolor <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if plot/x/y/width/height arguments are missing (could be in object)
  if(!hasArg(plot)) plot <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(width)) width <- NULL
  if(!hasArg(height)) height <- NULL

  ## Compile all parameters into an internal object
  bb_heatmapLegendInternal <- structure(list(plot = plot, border = border, x = x, y = y, width = width, height = height,
                                     orientation = orientation, fontsize = fontsize, fontcolor = fontcolor,
                                     just = just, default.units = default.units, gp = NULL), class = "bb_heatmapLegendInternal")

  bb_heatmapLegendInternal <- parseParams(bb_params = params, object_params = bb_heatmapLegendInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_heatmapLegendInternal$border)) bb_heatmapLegendInternal$border <- FALSE
  if(is.null(bb_heatmapLegendInternal$orientation)) bb_heatmapLegendInternal$orientation <- "v"
  if(is.null(bb_heatmapLegendInternal$fontsize)) bb_heatmapLegendInternal$fontsize <- 8
  if(is.null(bb_heatmapLegendInternal$fontcolor)) bb_heatmapLegendInternal$fontcolor <- "dark grey"
  if(is.null(bb_heatmapLegendInternal$just)) bb_heatmapLegendInternal$just <- c("left", "top")
  if(is.null(bb_heatmapLegendInternal$default.units)) bb_heatmapLegendInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # INITIAL ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_heatmapLegendInternal$plot)){
    stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  }

  if (is.null(bb_heatmapLegendInternal$plot$color_palette)){

    stop("Cannot add heatmap legend to an input plot that does not have a color palette.", call. = FALSE)

  }

  if (is.null(bb_heatmapLegendInternal$plot$zrange)){

    stop("Cannot add heatmap legend to an input plot that does not have a zrange.", call. = FALSE)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_heatmapLegend <- structure(list(color_palette = bb_heatmapLegendInternal$plot$color_palette, min_val = bb_heatmapLegendInternal$plot$zrange[1],
                                     max_val = bb_heatmapLegendInternal$plot$zrange[2], x = bb_heatmapLegendInternal$x, y = bb_heatmapLegendInternal$y, width = bb_heatmapLegendInternal$width,
                                     height = bb_heatmapLegendInternal$height, just = bb_heatmapLegendInternal$just, grobs = NULL), class = "bb_heatmapLegend")

  # ======================================================================================================================================================================================
  # CALL ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Must have a BentoBox page with a plot before adding a heatmap legend.")
  errorcheck_bb_annoHeatmapLegend(bb_heatmapLegend = bb_heatmapLegend, orientation = bb_heatmapLegendInternal$orientation)
  if(is.null(bb_heatmapLegend$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_heatmapLegend$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_heatmapLegend$width)) stop("argument \"width\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_heatmapLegend$height)) stop("argument \"height\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  bb_heatmapLegend <- defaultUnits(object = bb_heatmapLegend, default.units = bb_heatmapLegendInternal$default.units)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = bb_heatmapLegend)

  ## Make viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_heatmapLegend", length(grep(pattern = "bb_heatmapLegend", x = currentViewports)) + 1)

  ## Make viewport
  vp <- viewport(height = page_coords$height, width = page_coords$width,
                 x = page_coords$x, y = page_coords$y, just = bb_heatmapLegend$just, name = vp_name)

  pushViewport(vp)

  # ======================================================================================================================================================================================
  # SCALE COLORS
  # ======================================================================================================================================================================================

  color_scale <- bb_maptocolors(vec = seq(bb_heatmapLegend$min_val, bb_heatmapLegend$max_val, length.out = 100), col = bb_heatmapLegend$color_palette)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("heatmapLegend_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # CAPTURE AND SEPARATE "..." PARAMETERS
  # ======================================================================================================================================================================================

  params <- list(...)
  if (length(params) != 0){

    bb_heatmapLegendInternal$gp <- gpar(...)
    bb_heatmapLegendInternal$gp$fontsize <- bb_heatmapLegendInternal$fontsize
    bb_heatmapLegendInternal$gp$col <- bb_heatmapLegendInternal$fontcolor

    if ("col" %in% names(params)){
      bb_heatmapLegendInternal$gp$linecol <- params$col
    } else {
      bb_heatmapLegendInternal$gp$linecol <- "black"
    }

  } else {
    bb_heatmapLegendInternal$gp <- gpar(fontsize = bb_heatmapLegendInternal$fontsize, col = bb_heatmapLegendInternal$fontcolor)
  }

  # ======================================================================================================================================================================================
  # VERTICAL ORIENTATION
  # ======================================================================================================================================================================================
  if (bb_heatmapLegendInternal$orientation == "v"){

    digitLab <- textGrob(label = 0, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = bb_heatmapLegendInternal$gp)
    lowLab <- textGrob(label = bb_heatmapLegend$min_val, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = bb_heatmapLegendInternal$gp)
    highLab <- textGrob(label = bb_heatmapLegend$max_val, x = 0.5, y = 1, just = c("center", "top"), default.units = "npc",
                        gp = bb_heatmapLegendInternal$gp)

    lH <- convertHeight(x = grobHeight(lowLab), unitTo = "npc", valueOnly = T)
    hH <- convertHeight(x = grobHeight(highLab), unitTo = "npc", valueOnly = T)
    dH <- convertHeight(x = grobHeight(digitLab), unitTo = "npc", valueOnly = T)

    new_height <- 1 - (lH + hH + dH)

    color_scale <- rasterGrob(rev(color_scale), width = page_coords$width, height = unit(new_height, "npc"),
                              y = unit(1 - (hH + (0.5 * dH)), "npc"), x = unit(0.5, "npc"), just = "top")


    if (bb_heatmapLegendInternal$border == T){
      bb_heatmapLegendInternal$gp$fill <- NA
      bb_heatmapLegendInternal$gp$col <- bb_heatmapLegendInternal$gp$linecol
      borderGrob <- rectGrob(y = unit(1 - (hH + (0.5 * dH)), "npc"), just = "top", width = page_coords$width, height = unit(new_height, "npc"),
                gp = bb_heatmapLegendInternal$gp)


    }

  }

  # ======================================================================================================================================================================================
  # HORIZONTAL ORIENTATION
  # ======================================================================================================================================================================================
  if (bb_heatmapLegendInternal$orientation == "h"){

    digitLab <- textGrob(label = 0, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                         gp = bb_heatmapLegendInternal$gp)
    lowLab <- textGrob(label = bb_heatmapLegend$min_val, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                       gp = bb_heatmapLegendInternal$gp)
    highLab <- textGrob(label = bb_heatmapLegend$max_val, x = 1, y = 0.5, just = c("right", "center"), default.units = "npc",
                        gp = bb_heatmapLegendInternal$gp)

    lW <- convertWidth(x = grobWidth(lowLab), unitTo = "npc", valueOnly = T)
    hW <- convertWidth(x = grobWidth(highLab), unitTo = "npc", valueOnly = T)
    dW <- convertWidth(x = grobWidth(digitLab), unitTo = "npc", valueOnly = T)

    new_width <- 1 - (hW + lW + dW)

    color_scale <- rasterGrob(matrix(data = color_scale, nrow = 1, ncol = length(color_scale)), width = unit(new_width, "npc"), height = page_coords$height,
                x = unit(lW + (0.5 * dW), "npc"), just = "left")

    if (bb_heatmapLegendInternal$border == T){
      bb_heatmapLegendInternal$gp$fill <- NA
      bb_heatmapLegendInternal$gp$col <- bb_heatmapLegendInternal$gp$linecol
      borderGrob <- rectGrob(x = unit(lW + (0.5 * dW), "npc"), just = "left", width = unit(new_width, "npc"), height = page_coords$height,
                gp = bb_heatmapLegendInternal$gp)

    }

  }

  ## Go back to root viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ADD GROBS TO GTREE AND ASSIGN TO OBJECT
  # ======================================================================================================================================================================================

  ## Add grobs to gTree
  assign("heatmapLegend_grobs", setChildren(get("heatmapLegend_grobs", envir = bbEnv), children = gList(lowLab, highLab, color_scale)), envir = bbEnv)
  if (bb_heatmapLegendInternal$border == T){
    assign("heatmapLegend_grobs", addGrob(gTree = get("heatmapLegend_grobs", envir = bbEnv), child = borderGrob), envir = bbEnv)
  }

  ## Add grobs to object
  bb_heatmapLegend$grobs = get("heatmapLegend_grobs", envir = bbEnv)
  grid.draw(bb_heatmapLegend$grobs)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_heatmapLegend)

}
