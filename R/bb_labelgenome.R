#' @export

bb_labelgenome <- function(chrom, chromstart, chromend, scale = "bp", width = 3.25, height = 0.125, x = 0.75, y = 4, units = "inches",
                            fontsize = 10){

  if (scale == "bp"){
    fact = 1
  }

  if (scale == "Mb"){
    fact = 1000000
  }

  if (scale == "Kb"){
    fact = 1000
  }

  chromstartlabel = chromstart/fact
  chromendlabel = chromend/fact

  ##PLOT THE AXIS##

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)


  ## Convert x and y coordinates and height and width to same page_units
  old_x <- unit(x, units = units)
  old_y <- unit(y, units = units)
  old_height <- unit(height, units = units)
  old_width <- unit(width, units = units)
  new_x <- convertX(old_x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(old_y, unitTo = page_units, valueOnly = TRUE)
  new_height <- convertHeight(old_height, unitTo = page_units, valueOnly = TRUE)
  new_width <- convertWidth(old_width, unitTo = page_units, valueOnly = TRUE)

  ## Convert coordinates for viewport
  converted_coords = convert_coordinates(height = new_height, width = new_width, x = new_x, y = new_y, pageheight = page_height)

  vp <- viewport(width = unit(new_width, page_units), height = unit(new_height, page_units), x = unit(converted_coords[1], units = page_units),
                 y = unit(converted_coords[2], units = page_units))
  pushViewport(vp)

  grid.segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = "grey"))

  grid.text(chrom, x = 0.5, y = 0, gp = gpar(fontface = "bold", fontsize = fontsize, col = "grey"))
  grid.text(paste(chromstartlabel, scale, sep = " "), x = 0, y = 0, just = "left", gp = gpar(fontsize = fontsize, col = "grey"))
  grid.text(paste(chromendlabel, scale, sep = " "), x = 1.025, y = 0, just = "right", gp = gpar(fontsize = fontsize, col = "grey"))

  ## go back up a viewport
  upViewport()


}





