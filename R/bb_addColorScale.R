#' adds color scale for heatmap plots
#'
#' @param plot plot to add scale to
#' @param x A numeric vector or unit object specifying x-location
#' @param y A numeric vector or unit object specifying y-location
#' @param width A numeric vector or unit object specifying width
#' @param height A numeric vector or unit object specifying height
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param border option to add border around scale
#' @param orientation "v" (vertical) or "h" (horizontal) orientation
#' @param fontsize fontsize for text
#' @param fontcolor fontcolor for text
#' @param just justification of scale viewport
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#'
#' @author Nicole Kramer
#' @export
bb_addColorScale <- function(plot, x, y, width, height, params = NULL, border = FALSE, orientation = "v", fontsize = 8,
                      fontcolor = "dark grey", just = c("left", "top"), default.units = "inches", ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_addscale
  errorcheck_bb_addscale <- function(bb_scale){

    ## checking min_val and max val
    if (as.numeric(bb_scale$min_val) > as.numeric(bb_scale$max_val)){

      warning("/'min_val/' is larger than /'max_val/'. Scale labels may be incorrect.", call. = FALSE)

    }

    ## check proper orientation
    if (!bb_scale$orientation %in% c("v", "h")){

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
  bb_scaleInternal <- structure(list(plot = plot, border = border, x = x, y = y, width = width, height = height,
                                     orientation = orientation, fontsize = fontsize, fontcolor = fontcolor,
                                     just = just, default.units = default.units), class = "bb_scaleInternal")

  bb_scaleInternal <- parseParams(bb_params = params, object_params = bb_scaleInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_scaleInternal$border)) bb_scaleInternal$border <- FALSE
  if(is.null(bb_scaleInternal$orientation)) bb_scaleInternal$orientation <- "v"
  if(is.null(bb_scaleInternal$fontsize)) bb_scaleInternal$fontsize <- 8
  if(is.null(bb_scaleInternal$fontcolor)) bb_scaleInternal$fontcolor <- "dark grey"
  if(is.null(bb_scaleInternal$just)) bb_scaleInternal$just <- c("left", "top")
  if(is.null(bb_scaleInternal$default.units)) bb_scaleInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # INITIAL ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_scaleInternal$plot)){
    stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  }

  if (is.null(bb_scaleInternal$plot$color_palette)){

    stop("Cannot add color scale to an input plot that does not have a color palette.", call. = FALSE)

  }

  if (is.null(bb_scaleInternal$plot$zrange)){

    stop("Cannot add color scale to an input plot that does not have a zrange.", call. = FALSE)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_scale <- structure(list(color_palette = bb_scaleInternal$plot$color_palette, min_val = bb_scaleInternal$plot$zrange[1], max_val = bb_scaleInternal$plot$zrange[2],
                             orientation = bb_scaleInternal$orientation, height = bb_scaleInternal$height, width = bb_scaleInternal$width, x = bb_scaleInternal$x,
                             y = bb_scaleInternal$y, just = bb_scaleInternal$just, grobs = NULL,
                             gp = NULL), class = "bb_scale")

  # ======================================================================================================================================================================================
  # CALL ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Must have a BentoBox page with a plot before adding a scale.")
  errorcheck_bb_addscale(bb_scale = bb_scale)
  if(is.null(bb_scale$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_scale$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_scale$width)) stop("argument \"width\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_scale$height)) stop("argument \"height\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  bb_scale <- defaultUnits(object = bb_scale, default.units = bb_scaleInternal$default.units)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = bb_scale)

  ## Make viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_scale", length(grep(pattern = "bb_scale", x = currentViewports)) + 1)

  ## Make viewport
  vp <- viewport(height = page_coords$height, width = page_coords$width,
                 x = page_coords$x, y = page_coords$y, just = bb_scale$just, name = vp_name)

  pushViewport(vp)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("scale_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # CAPTURE AND SEPARATE "..." PARAMETERS
  # ======================================================================================================================================================================================

  params <- list(...)
  if (length(params) != 0){

    bb_scale$gp <- gpar(...)
    bb_scale$gp$fontsize <- bb_scaleInternal$fontsize
    bb_scale$gp$col <- bb_scaleInternal$fontcolor

    if ("col" %in% names(params)){
      bb_scale$gp$linecol <- params$col
    } else {
      bb_scale$gp$linecol <- "black"
    }

  } else {
    bb_scale$gp <- gpar(fontsize = bb_scaleInternal$fontsize, col = bb_scaleInternal$fontcolor)
  }

  # ======================================================================================================================================================================================
  # VERTICAL ORIENTATION
  # ======================================================================================================================================================================================
  if (bb_scaleInternal$orientation == "v"){

    digitLab <- textGrob(label = 0, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = bb_scale$gp)
    lowLab <- textGrob(label = bb_scale$min_val, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = bb_scale$gp)
    highLab <- textGrob(label = bb_scale$max_val, x = 0.5, y = 1, just = c("center", "top"), default.units = "npc",
                        gp = bb_scale$gp)

    lH <- convertHeight(x = grobHeight(lowLab), unitTo = "npc", valueOnly = T)
    hH <- convertHeight(x = grobHeight(highLab), unitTo = "npc", valueOnly = T)
    dH <- convertHeight(x = grobHeight(digitLab), unitTo = "npc", valueOnly = T)

    new_height <- 1 - (lH + hH + dH)

    color_scale <- rasterGrob(rev(bb_scale$color_palette), width = page_coords$width, height = unit(new_height, "npc"),
                              y = unit(1 - (hH + (0.5 * dH)), "npc"), x = unit(0.5, "npc"), just = "top")


    if (bb_scaleInternal$border == T){
      bb_scale$gp$fill <- NA
      bb_scale$gp$col <- bb_scale$gp$linecol
      borderGrob <- rectGrob(y = unit(1 - (hH + (0.5 * dH)), "npc"), just = "top", width = page_coords$width, height = unit(new_height, "npc"),
                gp = bb_scale$gp)


    }

  }

  # ======================================================================================================================================================================================
  # HORIZONTAL ORIENTATION
  # ======================================================================================================================================================================================
  if (bb_scaleInternal$orientation == "h"){

    digitLab <- textGrob(label = 0, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                         gp = bb_scale$gp)
    lowLab <- textGrob(label = bb_scale$min_val, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                       gp = bb_scale$gp)
    highLab <- textGrob(label = bb_scale$max_val, x = 1, y = 0.5, just = c("right", "center"), default.units = "npc",
                        gp = bb_scale$gp)

    lW <- convertWidth(x = grobWidth(lowLab), unitTo = "npc", valueOnly = T)
    hW <- convertWidth(x = grobWidth(highLab), unitTo = "npc", valueOnly = T)
    dW <- convertWidth(x = grobWidth(digitLab), unitTo = "npc", valueOnly = T)

    new_width <- 1 - (hW + lW + dW)

    color_scale <- rasterGrob(matrix(data = bb_scale$color_palette, nrow = 1, ncol = length(bb_scale$color_palette)), width = unit(new_width, "npc"), height = page_coords$height,
                x = unit(lW + (0.5 * dW), "npc"), just = "left")

    if (bb_scaleInternal$border == T){
      bb_scale$gp$fill <- NA
      bb_scale$gp$col <- bb_scale$gp$linecol
      borderGrob <- rectGrob(x = unit(lW + (0.5 * dW), "npc"), just = "left", width = unit(new_width, "npc"), height = page_coords$height,
                gp = bb_scale$gp)

    }

  }

  ## Go back to root viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ADD GROBS TO GTREE AND ASSIGN TO SCALE OBJECT
  # ======================================================================================================================================================================================

  ## Add grobs to gTree
  assign("scale_grobs", setChildren(get("scale_grobs", envir = bbEnv), children = gList(lowLab, highLab, color_scale)), envir = bbEnv)
  if (bb_scaleInternal$border == T){
    assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = borderGrob), envir = bbEnv)
  }

  ## Add grobs to scale object
  bb_scale$grobs = get("scale_grobs", envir = bbEnv)
  grid.draw(bb_scale$grobs)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_scale)

}
