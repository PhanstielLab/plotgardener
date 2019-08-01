#' adds legend to HiC interaction matrix
#'
#' @param color_vector sequential list of colors used to represent interaction frequency; returned from bb_hic function
#' @param min_label label for minimum value, i.e first value of zrange
#' @param max_label label for maximum value, i.e. last value of zrange
#' @param height legend height, in inches
#' @param width legend width, in inches
#' @param x legend x-coordinate
#' @param y legend y-coordinate
#' @param pageheight height of page on which plot will be placed
#' @param orientation "vertical" or "horizontal" orientation
#' @param fontsize fontsize for text
#' @param fontcolor fontcolor for test
#'
#' @author Nicole Kramer
#' @export

bb_legend <- function(color_vector, min_label, max_label, height, width, x, y, pageheight, orientation = "vertical", fontsize = 8, fontcolor = "grey", ...){
  ## Get page_height and margins from bbEnv
  page_height <- get("height", envir = bbEnv)
  top_margin <- get("top_margin", envir = bbEnv)
  left_margin <- get("left_margin", envir = bbEnv)

  x = x + left_margin
  y = y + top_margin

  ## Convert user-inputted coordinates
  legend_coordinates <- convert_coordinates(height = height, width = width, x = x, y = y, pageheight = page_height)

  # VIEWPORTS
  # ======================================================================================================================================================================================
  ## Go up a viewport if not at the root of all viewports
  if(is.null(current.vpPath()) == FALSE){
    upViewport()
  }

  vp <- viewport(height = unit(height, "in"), width = unit(width, "in"), x = unit(legend_coordinates[1], "in"), y = unit(legend_coordinates[2], "in"))
  pushViewport(vp)
  # ======================================================================================================================================================================================

  ## Get length of color vector to determine how to divide the space
  ncolors <- length(color_vector)

  # VERTICAL ORIENTATION
  # ======================================================================================================================================================================================
  if (orientation == "vertical"){

    ## All rectangles have same x-coordinate, and increase in y-coordinate by the height of each rectangle
    xcoords <- rep(0.5, ncolors)
    ycoords <- seq(from = (1/ncolors)/2, to = 1, by = 1/ncolors)
    grid.rect(xcoords, ycoords, height = 1/ncolors, width = 1, gp = gpar(col = color_vector, fill = color_vector))

    ## Min label goes at the bottom, max label goes at the top
    grid.text(min_label, x = 0.5, y = -0.1, gp = gpar(fontsize = fontsize, col = fontcolor))
    grid.text(max_label, x = 0.5, y = 1.1, gp = gpar(fontsize = fontsize, col = fontcolor))
  }

  # HORIZONTAL ORIENTATION
  # ======================================================================================================================================================================================
  if (orientation == "horizontal"){

    ## All rectanges have same y-coordinate, and increase in x-coordinate by the width of each rectangle
    xcoords <- seq(from = (1/ncolors)/2, to = 1, by = 1/ncolors)
    ycoords <- rep(0.5, ncolors)
    grid.rect(xcoords, ycoords, height = 1, width = 1/ncolors, gp = gpar(col = color_vector, fill = color_vector))

    ## Min label goes on left, max label goes on right
    grid.text(min_label, x = -0.1, y = 0.5, gp = gpar(fontsize = fontsize, col = fontcolor))
    grid.text(max_label, x = 1.1, y = 0.5, gp = gpar(fontsize = fontsize, col = fontcolor))
  }

  ## Go back to root viewport
  upViewport()
}
