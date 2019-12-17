#' adds color scale for heatmap plots
#'
#' @param plot plot to add scale to
#' @param border option to add border around scale
#' @param height A unit object specifying height
#' @param width A unit object specifying width
#' @param x A unit object specifying x-location
#' @param y A unit object specifying y-location
#' @param orientation "vertical" or "horizontal" orientation
#' @param fontsize fontsize for text
#' @param fontcolor fontcolor for text
#' @param fontface fontface for text
#' @param just justification of scale viewport
#'
#' @author Nicole Kramer
#' @export

bb_addScale <- function(plot, border = FALSE, height, width, x, y, orientation = "vertical", fontsize = 8,
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

  bb_scale <- structure(list(color_palette = plot$color_palette, min_val = plot$zrange[1], max_val = plot$zrange[2],
                              orientation = orientation, height = height, width = width, x = x, y = y, grobs = NULL,
                             gpar = list(fontsize = fontsize, fontcolor = fontcolor, fontface = fontface)), class = "bb_scale")

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
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_scale", length(grep(pattern = "bb_scale", x = current_viewports)) + 1)

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
  if (orientation == "vertical"){

    digitLab <- textGrob(label = 0, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    lowLab <- grid.text(label = min_val, x = 0.5, y = 0, just = c("center", "bottom"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    highLab <- grid.text(label = max_val, x = 0.5, y = 1, just = c("center", "top"), default.units = "npc",
                        gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))

    lH <- convertHeight(x = grobHeight(lowLab), unitTo = "npc", valueOnly = T)
    hH <- convertHeight(x = grobHeight(highLab), unitTo = "npc", valueOnly = T)
    dH <- convertHeight(x = grobHeight(digitLab), unitTo = "npc", valueOnly = T)

    new_height <- 1 - (lH + hH + dH)

    color_scale <- rasterGrob(rev(color_vector), width = unit(page_coords[[1]]$width, page_coords[[3]]), height = unit(new_height, "npc"), y = 1 - (hH + (0.5 * dH)), just = "top")
    grid.draw(color_scale)

    if (border == T){

      borderGrob <- grid.rect(y = 1 - (hH + (0.5 * dH)), just = "top", width = unit(page_coords[[1]]$width, page_coords[[3]]), height = unit(new_height, "npc"),
                gp = gpar(col = "black", fill = NA) )


    }

  }

  # ======================================================================================================================================================================================
  # HORIZONTAL ORIENTATION
  # ======================================================================================================================================================================================
  if (orientation == "horizontal"){

    digitLab <- textGrob(label = 0, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                         gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    lowLab <- grid.text(label = min_val, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc",
                       gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    highLab <- grid.text(label = max_val, x = 1, y = 0.5, just = c("right", "center"), default.units = "npc",
                        gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))

    lW <- convertWidth(x = grobWidth(lowLab), unitTo = "npc", valueOnly = T)
    hW <- convertWidth(x = grobWidth(highLab), unitTo = "npc", valueOnly = T)
    dW <- convertWidth(x = grobWidth(digitLab), unitTo = "npc", valueOnly = T)

    new_width <- 1 - (hW + lW + dW)

    color_scale <- rasterGrob(matrix(data = color_vector, nrow = 1, ncol = length(color_vector)), width = unit(new_width, "npc"), height = unit(page_coords[[1]]$height, page_coords[[3]]),
                x = lW + (0.5 * dW), just = "left")
    grid.draw(color_scale)

    if (border == T){

      borderGrob <- grid.rect(x = lW + (0.5 * dW), just = "left", width = unit(new_width, "npc"), height = unit(page_coords[[1]]$height, page_coords[[3]]),
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

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_scale)

}
