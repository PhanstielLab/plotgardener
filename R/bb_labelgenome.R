#' Adds genome coordinates to the axis of a plot
#'
#' @param plot plot to annotate (saved as output from BentoBox plotting function)
#' @param scale scale of the plot; options are "bp", "Kb", or "Mb"
#' @param side side of the plot to annotate, options are "bottom", "top", "right", or "left"
#' @param fontsize fontsize for labels (in points)
#' @param color color for all text/lines


#' @export

bb_labelgenome <- function(plot, scale = "bp", side = "bottom", fontsize = 10, color = "grey"){


  #chrom
  #chromstart
  #chromend
  #width
  #height
  #x
  #y
  #units

  ## Grab information from input plot
  chrom <- plot$chrom
  chromstart <- plot$chromstart
  chromend <- plot$chromend
  plot_height <- plot$height
  plot_width <- plot$width
  plot_x <- plot$x
  plot_y <- plot$y
  plot_units <- plot$units


  ## Determine scale of labels
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

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  ## Make and place viewport based on input plot

  if (side == "bottom"){

    label_width <- convertWidth(unit(plot_width, units = plot_units), unitTo = page_units, valueOnly = TRUE)
    label_height <- convertHeight(unit(fontsize/72, units = "inches"), unitTo = page_units, valueOnly = TRUE)

    label_y <- convertY(unit((plot_y + plot_height), units = plot_units), unitTo = page_units, valueOnly = TRUE)
    label_x <- convertX(unit(plot_x, units = plot_units), unitTo = page_units, valueOnly = TRUE)

    converted_coords = convert_coordinates(height = label_height, width = label_width, x = label_x, y = label_y, pageheight = page_height)
    vp <- viewport(width = unit(label_width, page_units), height = unit(label_height, page_units), x = unit(converted_coords[1], units = page_units),
                   y = unit(converted_coords[2], units = page_units))
    pushViewport(vp)

    grid.segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = color))
    grid.text(chrom, x = 0.5, y = 0.25, gp = gpar(fontface = "bold", fontsize = fontsize, col = color))
    grid.text(paste(chromstartlabel, scale, sep = " "), x = 0, y = 0.25, just = "left", gp = gpar(fontsize = fontsize, col = color))
    grid.text(paste(chromendlabel, scale, sep = " "), x = 1, y = 0.25, just = "right", gp = gpar(fontsize = fontsize, col = color))

  } else if (side == "top"){

    label_width <- convertWidth(unit(plot_width, units = plot_units), unitTo = page_units, valueOnly = TRUE)
    label_height <- convertHeight(unit(fontsize/72, units = "inches"), unitTo = page_units, valueOnly = TRUE)
    label_height2 <- convertHeight(unit(label_height, units = page_units), unitTo = plot_units, valueOnly = TRUE)

    label_y <- convertY(unit((plot_y - label_height2), units = plot_units), unitTo = page_units, valueOnly = TRUE)
    label_x <- convertX(unit(plot_x, units = plot_units), unitTo = page_units, valueOnly = TRUE)

    converted_coords = convert_coordinates(height = label_height, width = label_width, x = label_x, y = label_y, pageheight = page_height)
    vp <- viewport(width = unit(label_width, page_units), height = unit(label_height, page_units), x = unit(converted_coords[1], units = page_units),
                   y = unit(converted_coords[2], units = page_units))
    pushViewport(vp)

    grid.segments(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(col = color))
    grid.text(chrom, x = 0.5, y = 0.75, gp = gpar(fontface = "bold", fontsize = fontsize, col = color))
    grid.text(paste(chromstartlabel, scale, sep = " "), x = 0, y = 0.75, just = "left", gp = gpar(fontsize = fontsize, col = color))
    grid.text(paste(chromendlabel, scale, sep = " "), x = 1, y = 0.75, just = "right", gp = gpar(fontsize = fontsize, col = color))

  } else if (side == "left"){

    label_height <- convertHeight(unit(plot_width, units = plot_units), unitTo = page_units, valueOnly = TRUE)
    label_width <- convertWidth(unit(fontsize/72, units = "inches"), unitTo = page_units, valueOnly = TRUE)
    label_width2 <- convertWidth(unit(label_width, units = page_units), unitTo = plot_units, valueOnly = TRUE)

    label_y <- convertY(unit(plot_y, units = plot_units), unitTo = page_units, valueOnly = TRUE)
    label_x <- convertX(unit((plot_x - label_width2), units = plot_units), unitTo = page_units, valueOnly = TRUE)

    converted_coords = convert_coordinates(height = label_height, width = label_width, x = label_x, y = label_y, pageheight = page_height)
    vp <- viewport(width = unit(label_width, page_units), height = unit(label_height, page_units), x = unit(converted_coords[1], units = page_units),
                   y = unit(converted_coords[2], units = page_units))
    pushViewport(vp)

    grid.segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, gp = gpar(col = color))
    grid.text(chrom, x = 0.25, y = 0.5, gp = gpar(fontface = "bold", fontsize = fontsize, col = color), rot = 90)
    grid.text(paste(chromstartlabel, scale, sep = " "), x = 0.25, y = 0, just = "left", gp = gpar(fontsize = fontsize, col = color), rot = 90)
    grid.text(paste(chromendlabel, scale, sep = " "), x = 0.25, y = 1, just = "right", gp = gpar(fontsize = fontsize, col = color), rot = 90)

  } else if (side == "right"){

    label_height <- convertHeight(unit(plot_width, units = plot_units), unitTo = page_units, valueOnly = TRUE)
    label_width <- convertWidth(unit(fontsize/72, units = "inches"), unitTo = page_units, valueOnly = TRUE)

    label_y <- convertY(unit(plot_y , units = plot_units), unitTo = page_units, valueOnly = TRUE)
    label_x <- convertX(unit((plot_x + plot_width), units = plot_units), unitTo = page_units, valueOnly = TRUE)

    converted_coords = convert_coordinates(height = label_height, width = label_width, x = label_x, y = label_y, pageheight = page_height)
    vp <- viewport(width = unit(label_width, page_units), height = unit(label_height, page_units), x = unit(converted_coords[1], units = page_units),
                   y = unit(converted_coords[2], units = page_units))
    pushViewport(vp)

    grid.segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, gp = gpar(col = color))
    grid.text(chrom, x = 0.75, y = 0.5, gp = gpar(fontface = "bold", fontsize = fontsize, col = color), rot = 90)
    grid.text(paste(chromstartlabel, scale, sep = " "), x = 0.75, y = 0, just = "left", gp = gpar(fontsize = fontsize, col = color), rot = 90)
    grid.text(paste(chromendlabel, scale, sep = " "), x = 0.75, y = 1, just = "right", gp = gpar(fontsize = fontsize, col = color), rot = 90)

  }

  ## Go back up a viewport
  upViewport()


}





