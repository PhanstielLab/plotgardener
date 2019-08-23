#' adds legend to HiC interaction matrix
#'
#' @param color_vector sequential list of colors used to represent interaction frequency; returned from bb_hic function
#' @param min_label label for minimum value, i.e first value of zrange
#' @param max_label label for maximum value, i.e. last value of zrange
#' @param height legend height, in inches
#' @param width legend width, in inches
#' @param x legend x-coordinate
#' @param y legend y-coordinate
#' @param units units of height, width, and x and y coordinates
#' @param orientation "vertical" or "horizontal" orientation
#' @param fontsize fontsize for text
#' @param fontcolor fontcolor for test
#'
#' @author Nicole Kramer
#' @export

bb_legend <- function(color_vector, min_label, max_label, height, width, x, y, units = "inches", orientation = "vertical", fontsize = 8, fontcolor = "grey", ...){

  ## Get page_height from bbEnv
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  ## Convert x and y coordinates and height and width to same page_units
  old_x <- unit(x, units = units)
  old_y <- unit(y, units = units)
  old_width <- unit(width, units = units)
  old_height <- unit(height, units = units)
  new_x <- convertX(old_x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(old_y, unitTo = page_units, valueOnly = TRUE)
  new_height <- convertHeight(old_height, unitTo = page_units, valueOnly = TRUE)
  new_width <- convertWidth(old_width, unitTo = page_units, valueOnly = TRUE)

  ## Convert user-inputted coordinates
  legend_coordinates <- convert_coordinates(height = new_height, width = new_width, x = new_x, y = new_y, pageheight = page_height)

  # VIEWPORTS
  # ======================================================================================================================================================================================

  vp <- viewport(height = unit(height, page_units), width = unit(width, page_units), x = unit(legend_coordinates[1], page_units), y = unit(legend_coordinates[2], page_units))
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
