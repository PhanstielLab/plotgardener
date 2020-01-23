#' wrapper to draw a textGrob based on bb_makePage coordinates and units
#'
#' @param label character or expression
#' @param x A unit object specifying x-location
#' @param y A unit object specifying y-location
#' @param just justification of text relative to its (x, y) location
#' @param rotation angle to rotate text
#' @param fontcolor fontcolor
#' @param transparency degree of text transparency
#' @param fontsize the size of text (in points)
#' @param cex multiplier applied to fontsize
#' @param fontfamily the font family
#' @param fontface the fontface (bold, italic, ...)
#'
#' @export
bb_addText <- function(label, x, y, just = "center", rotation = 0, fontcolor = "black", transparency = 1, fontsize = 12, cex = 1,
                    fontfamily = "", fontface = "plain"){

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot add text without a BentoBox page.")

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_text <- structure(list(label = label, x = x, y = y, rotation = rotation, just = just, grobs = NULL,
                            gpar = list(fontcolor = fontcolor, transparency = transparency, fontsize = fontsize, cex = cex, fontfamily = fontfamily,
                                        fontface = fontface)), class = "bb_text")

  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  ## Convert coordinates to page_units
  new_x <- convertX(x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(y, unitTo = page_units, valueOnly = TRUE)

  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  text <- grid.text(label = label, x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), just = just, rot = rotation,
            gp = gpar(col = fontcolor, alpha = transparency, fontsize = fontsize, cex = cex, fontfamily = fontfamily, fontface = fontface))

  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================

  bb_text$grobs <- text

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_text)
}
