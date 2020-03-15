#' wrapper to draw a textGrob based on bb_makePage coordinates and units
#'
#' @param label character or expression
#' @param x A numeric vector or unit object specifying x-location
#' @param y A numeric vector or unit object specifying y-location
#' @param just justification of text relative to its (x, y) location
#' @param rotation angle to rotate text
#' @param fontcolor fontcolor
#' @param transparency degree of text transparency
#' @param fontsize the size of text (in points)
#' @param cex multiplier applied to fontsize
#' @param fontfamily the font family
#' @param fontface the fontface (bold, italic, ...)
#' @param default.units A string indicating the default units to use if x or y are only given as numeric vectors
#'
#' @export
bb_addText <- function(label, x, y, just = "center", rotation = 0, fontcolor = "black", transparency = 1, fontsize = 12, cex = 1,
                    fontfamily = "", fontface = "plain", default.units = "inches"){

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

  if (class(x) != "unit"){

    if (!is.numeric(x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_text$x <- unit(x, default.units)

  }

  if (class(y) != "unit"){

    if (!is.numeric(y)){

      stop("y-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(default.units)){

      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_text$y <- unit(y, default.units)

  }

  ## Convert coordinates to page_units
  new_x <- convertX(bb_text$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(bb_text$y, unitTo = page_units, valueOnly = TRUE)

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
