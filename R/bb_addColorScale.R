#' adds color scale for heatmap plots
#'
#' @param plot plot to add scale to
#' @param border option to add border around scale
#' @param height A unit object specifying height
#' @param width A unit object specifying width
#' @param x A unit object specifying x-location
#' @param y A unit object specifying y-location
#' @param orientation "v" (vertical) or "h" (horizontal) orientation
#' @param fontsize fontsize for text
#' @param fontcolor fontcolor for text
#' @param fontfamily fontfamily for text
#' @param fontface fontface for text
#' @param just justification of scale viewport
#'
#' @author Nicole Kramer
#' @export
bb_addColorScale <- function(plot, border = FALSE, height, width, x, y, orientation = "v", fontsize = 8,
                      fontcolor = "dark grey", fontfamily = "", fontface = "plain", just = c("left", "top"), ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_addlegend
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
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_scale <- structure(list(color_palette = plot$color_palette, min_val = plot$zrange[1], max_val = plot$zrange[2],
                              orientation = orientation, height = height, width = width, x = x, y = y, grobs = NULL,
                             gpar = list(fontsize = fontsize, fontcolor = fontcolor, fontface = fontface, fontfamily = fontfamily)), class = "bb_scale")

  # ======================================================================================================================================================================================
  # CALL ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Must have a BentoBox page with a plot before adding a scale.")
  errorcheck_bb_addscale(bb_scale = bb_scale)

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
                 x = page_coords$x, y = page_coords$y, just = just, name = vp_name)

  pushViewport(vp)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("scale_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # VERTICAL ORIENTATION
  # ======================================================================================================================================================================================
  if (orientation == "v"){

    digitLab <- textGrob(label = 0, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface, fontfamily = fontfamily))
    lowLab <- textGrob(label = bb_scale$min_val, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface, fontfamily = fontfamily))
    highLab <- textGrob(label = bb_scale$max_val, x = 0.5, y = 1, just = c("center", "top"), default.units = "npc",
                        gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface, fontfamily = fontfamily))

    lH <- convertHeight(x = grobHeight(lowLab), unitTo = "npc", valueOnly = T)
    hH <- convertHeight(x = grobHeight(highLab), unitTo = "npc", valueOnly = T)
    dH <- convertHeight(x = grobHeight(digitLab), unitTo = "npc", valueOnly = T)

    new_height <- 1 - (lH + hH + dH)

    color_scale <- rasterGrob(rev(bb_scale$color_palette), width = page_coords$width, height = unit(new_height, "npc"),
                              y = unit(1 - (hH + (0.5 * dH)), "npc"), x = unit(0.5, "npc"), just = "top")


    if (border == T){

      borderGrob <- rectGrob(y = unit(1 - (hH + (0.5 * dH)), "npc"), just = "top", width = page_coords$width, height = unit(new_height, "npc"),
                gp = gpar(col = "black", fill = NA) )


    }

  }

  # ======================================================================================================================================================================================
  # HORIZONTAL ORIENTATION
  # ======================================================================================================================================================================================
  if (orientation == "h"){

    digitLab <- textGrob(label = 0, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                         gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface, fontfamily = fontfamily))
    lowLab <- textGrob(label = bb_scale$min_val, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface, fontfamily = fontfamily))
    highLab <- textGrob(label = bb_scale$max_val, x = 1, y = 0.5, just = c("right", "center"), default.units = "npc",
                        gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface, fontfamily = fontfamily))

    lW <- convertWidth(x = grobWidth(lowLab), unitTo = "npc", valueOnly = T)
    hW <- convertWidth(x = grobWidth(highLab), unitTo = "npc", valueOnly = T)
    dW <- convertWidth(x = grobWidth(digitLab), unitTo = "npc", valueOnly = T)

    new_width <- 1 - (hW + lW + dW)

    color_scale <- rasterGrob(matrix(data = bb_scale$color_palette, nrow = 1, ncol = length(bb_scale$color_palette)), width = unit(new_width, "npc"), height = page_coords$height,
                x = unit(lW + (0.5 * dW), "npc"), just = "left")

    if (border == T){

      borderGrob <- rectGrob(x = unit(lW + (0.5 * dW), "npc"), just = "left", width = unit(new_width, "npc"), height = page_coords$height,
                gp = gpar(col = "black", fill = NA) )

    }

  }

  ## Go back to root viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ADD GROBS TO GTREE AND ASSIGN TO SCALE OBJECT
  # ======================================================================================================================================================================================

  ## Add grobs to gTree
  assign("scale_grobs", setChildren(get("scale_grobs", envir = bbEnv), children = gList(lowLab, highLab, color_scale)), envir = bbEnv)
  if (border == T){
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
