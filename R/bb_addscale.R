#' adds color scale for heatmap plots
#'
#' @param color_vector sequential list of colors used to represent interaction frequency; can be returned from bb_plothic function
#' @param min_val label for minimum value, i.e first value of zrange
#' @param max_val label for maximum value, i.e. last value of zrange
#' @param border option to add border around scale
#' @param height scale height
#' @param width scale width
#' @param x scale x-coordinate
#' @param y scale y-coordinate
#' @param units units of height, width, and x and y coordinates
#' @param orientation "vertical" or "horizontal" orientation
#' @param fontsize fontsize for text
#' @param fontcolor fontcolor for text
#' @param fontface fontface for text
#' @param just justification of scale viewport
#'
#' @author Nicole Kramer
#' @export

bb_addscale <- function(color_vector, min_val, max_val, border = FALSE, height, width, x, y, units = "inches", orientation = "vertical", fontsize = 8,
                      fontcolor = "grey", fontface = "plain", just = c("left", "top"), ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_addlegend
  errorcheck_bb_addscale <- function(bb_scale){

    ## checking min_val and max val
    if (as.numeric(bb_scale$min_val) > as.numeric(bb_scale$max_val)){

      warning("/'min_val/' is larger than /'max_val/'. Scale labels may be incorrect.")

    }

    ## check proper orientation
    if (!bb_scale$orientation %in% c("vertical", "horizontal")){

      stop("Invalid /'orientation/' parameter. Options are 'vertical' or 'horizontal'.")

    }

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_scale <- structure(list(color_palette = color_vector, min_val = min_val, max_val = max_val,
                              orientation = orientation, height = height, width = width, x = x, y = y, units = units, grobs = NULL), class = "bb_scale")

  # ======================================================================================================================================================================================
  # CALL ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_addscale(bb_scale = bb_scale)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = bb_scale)

  ## Make viewport
  vp <- viewport(height = unit(page_coords[[1]]$height, page_coords[[3]]), width = unit(page_coords[[1]]$width, page_coords[[3]]),
                 x = unit(page_coords[[1]]$x, page_coords[[3]]), y = unit((page_coords[[2]]-page_coords[[1]]$y), page_coords[[3]]), just = just, name = "bb_scale")

  pushViewport(vp)

  # ======================================================================================================================================================================================
  # GROBS
  # ======================================================================================================================================================================================

  assign("scale_grobs", gTree(name = "scale_grobs"), envir = bbEnv)

  # ======================================================================================================================================================================================
  # VERTICAL ORIENTATION
  # ======================================================================================================================================================================================
  if (orientation == "vertical"){

    digitLab <- textGrob(label = 0, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    lowLab <- textGrob(label = min_val, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    highLab <- textGrob(label = max_val, x = 0.5, y = 1, just = c("center", "top"), default.units = "npc",
                        gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))

    lH <- convertHeight(x = grobHeight(lowLab), unitTo = "npc", valueOnly = T)
    hH <- convertHeight(x = grobHeight(highLab), unitTo = "npc", valueOnly = T)
    dH <- convertHeight(x = grobHeight(digitLab), unitTo = "npc", valueOnly = T)

    new_height <- 1 - (lH + hH + dH)

    color_scale <- rasterGrob(rev(color_vector), width = unit(page_coords[[1]]$width, page_coords[[3]]), height = unit(new_height, "npc"), y = 1 - (hH + (0.5 * dH)), just = "top")

    grid.draw(lowLab)
    grid.draw(highLab)
    grid.draw(color_scale)

    assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = lowLab), envir = bbEnv)
    assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = highLab), envir = bbEnv)
    assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = color_scale), envir = bbEnv)

    if (border == T){

      borderGrob <- rectGrob(y = 1 - (hH + (0.5 * dH)), just = "top", width = unit(page_coords[[1]]$width, page_coords[[3]]), height = unit(new_height, "npc"),
                gp = gpar(col = "black", fill = NA) )
      grid.draw(borderGrob)
      assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = borderGrob), envir = bbEnv)

    }

  }

  # ======================================================================================================================================================================================
  # HORIZONTAL ORIENTATION
  # ======================================================================================================================================================================================
  if (orientation == "horizontal"){

    digitLab <- textGrob(label = 0, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                         gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    lowLab <- textGrob(label = min_val, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    highLab <- textGrob(label = max_val, x = 1, y = 0.5, just = c("right", "center"), default.units = "npc",
                        gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))

    lW <- convertWidth(x = grobWidth(lowLab), unitTo = "npc", valueOnly = T)
    hW <- convertWidth(x = grobWidth(highLab), unitTo = "npc", valueOnly = T)
    dW <- convertWidth(x = grobWidth(digitLab), unitTo = "npc", valueOnly = T)

    new_width <- 1 - (hW + lW + dW)

    color_scale <- rasterGrob(matrix(data = color_vector, nrow = 1, ncol = length(color_vector)), width = unit(new_width, "npc"), height = unit(page_coords[[1]]$height, page_coords[[3]]),
                x = lW + (0.5 * dW), just = "left")

    grid.draw(lowLab)
    grid.draw(highLab)
    grid.draw(color_scale)

    assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = lowLab), envir = bbEnv)
    assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = highLab), envir = bbEnv)
    assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = color_scale), envir = bbEnv)

    if (border == T){

      borderGrob <- rectGrob(x = lW + (0.5 * dW), just = "left", width = unit(new_width, "npc"), height = unit(page_coords[[1]]$height, page_coords[[3]]),
                gp = gpar(col = "black", fill = NA) )
      grid.draw(borderGrob)
      assign("scale_grobs", addGrob(gTree = get("scale_grobs", envir = bbEnv), child = borderGrob), envir = bbEnv)

    }

  }

  grob_list <- lapply(get("scale_grobs", envir = bbEnv)$children, convert_gpath)
  bb_scale$grobs = grob_list

  ## Go back to root viewport
  upViewport()

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_scale)

}
