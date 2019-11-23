#' Adds genome coordinates to the axis of a plot
#'
#' @param plot plot to annotate (saved as output from BentoBox plotting function)
#' @param scale scale of the plot; options are "bp", "Kb", or "Mb"
#' @param side side of the plot to annotate, options are "bottom", "top", "right", or "left"
#' @param fontsize fontsize for labels (in points)
#' @param color color for all text/lines


#' @export

bb_labelgenome <- function(plot = NULL, chrom, chromstart, chromend, width, x, y, just = c("left", "top"), rotation = 0,
                           units = "inches", scale = "bp", fontsize = 10, fontcolor = "black", lineAbove = T, linecolor = "black", lwd = 1){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_labelgenome
  # errorcheck_bb_labelgenome <- function()

  ## Define a function to convert plot x and y into center of plot based on justification
  adjust_coords <- function(plot, page_units, page_height){

    plot_y <- convertY(unit(plot$y, units = plot$units), unitTo = page_units, valueOnly = TRUE)
    plot_y <- page_height - plot_y
    plot_y <- convertY(unit(plot_y, units = page_units), unitTo = plot$units, valueOnly = TRUE)

    if (length(plot$just == 2)){

      if ("left" %in% plot$just & "center" %in% plot$just){
        ## convert the x-coordinate only
        plot_x <- plot$x + (0.5 * plot$width)
      } else if ("right" %in% plot$just & "center" %in% plot$just){
        ## convert the x-coordinate only
        plot_x <- plot$x - (0.5 * plot$width)
      } else if ("center" %in% plot$just & "bottom" %in% plot$just){
        ## convert the y-coordinate only
        plot_x <- plot$x
        plot_y <- plot_y + (0.5 * plot$height)
      } else if ("center" %in% plot$just & "top" %in% plot$just){
        ## convert the y-coordinate only
        plot_x <- plot$x
        plot_y <- plot_y - (0.5 * plot$height)
      } else if ("left" %in% plot$just & "top" %in% plot$just){
        ## convert x-coordinate and y-coordinate
        plot_x <- plot$x + (0.5 * plot$width)
        plot_y <- plot_y - (0.5 * plot$height)
      } else if ("right" %in% plot$just & "top" %in% plot$just){
        ## convert x-coordinate and y-coordinate
        plot_x <- plot$x - (0.5 * plot$width)
        plot_y <- plot_y - (0.5 * plot$height)
      } else if ("left" %in% plot$just & "bottom" %in% plot$just){
        ## convert x-coordinate and y-coordinate
        plot_x <- plot$x + (0.5 * plot$width)
        plot_y <- plot_y + (0.5 * plot$height)
      } else if ("right" %in% plot$just & "bottom" %in% plot$just){
        ## convert x-coordinate and y-coordinate
        plot_x <- plot$x - (0.5 * plot$width)
        plot_y <- plot_y + (0.5 * plot$height)
      } else {
        ## no conversion
        plot_x <- plot$x
      }

    } else if (length(plot$just == 1)){

      if (plot$just == "left"){
        ## convert the x-coordinate only
        plot_x <- plot$x + (0.5 * plot$width)
      } else if (plot$just == "right"){
        ## convert the x-coordinate only
        plot_x <- plot$x - (0.5 * plot$width)
      } else if (plot$just == "bottom"){
        ## convert the y-coordinate only
        plot_x <- plot$x
        plot_y <- plot_y + (0.5 * plot$height)
      } else if (plot$just == "top"){
        ## convert the y-coordinate only
        plot_x <- plot$x
        plot_y <- plot_y - (0.5 * plot$height)
      } else {
        ## no conversion
        plot_x <- plot$x
      }

    }

    return(list(plot_x, plot_y))

  }

  ## Define a function to get label x and y based on plot input coordinates
  label_coords <- function(plot, xCoord, yCoord, side, genome_label){

    if (side == "bottom"){

      label_x <- xCoord - (0.5 * plot$width)
      label_y <- yCoord - (0.5 * plot$height)
      justification <- c("left", "top")

    } else if (side == "top"){

      label_x <- xCoord - (0.5 * plot$width)
      label_y <- yCoord + (0.5 * plot$height)
      justification <- c("left", "bottom")

    } else if (side == "left"){

      label_x <- xCoord - (0.5 * plot$width)
      label_y <- yCoord + (0.5 * plot$height)
      justification <- c("right", "top")

    } else if (side == "right"){

      label_x <- xCoord + (0.5 * plot$width)
      label_y <- yCoord + (0.5 * plot$height)
      justification <- c("left", "top")

    }

    genome_label[["x"]] = label_x
    genome_label[["y"]] = label_y
    genome_label[["units"]] = plot$units
    genome_label[["just"]] = justification

    return(genome_label)

  }

  ## Define a function to get label width and height based on side
  label_dims <- function(plot, side, text, page_units, genome_label){

    if (side == "bottom" | side == "top"){

      label_width <- plot$width
      label_height <- convertHeight(heightDetails(text), unitTo = plot$units, valueOnly = TRUE)

    } else if (side == "left" | side == "right"){

      label_width <- convertWidth(heightDetails(text), unitTo = plot$units, valueOnly = TRUE)
      label_height <- plot$height

    }

    genome_label[["width"]] = label_width
    genome_label[["height"]] = label_height

    return(genome_label)

  }

  ## Define a function to plot automatic label
  plot_label <- function(side, genome_label){

    if (side == "bottom"){

      grid.segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = genome_label$linecolor, lwd = genome_label$lwd))
      grid.text(label = genome_label$chromlabel, x = 0.5, y = 0.85, gp = gpar(fontface = "bold", fontsize = genome_label$fontsize, col = genome_label$fontcolor), just = c("center", "top"))
      grid.text(label = paste(round(genome_label$chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.85, just = c("left", "top"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor))
      grid.text(label = paste(round(genome_label$chromendlabel, 1), scale, sep = " "), x = 1, y = 0.85, just = c("right","top"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor))

    } else if (side == "top"){

      grid.segments(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(col = genome_label$linecolor, lwd = genome_label$lwd))
      grid.text(label = genome_label$chromlabel, x = 0.5, y = 0.15, gp = gpar(fontface = "bold", fontsize = genome_label$fontsize, col = genome_label$fontcolor), just = c("center", "bottom"))
      grid.text(label = paste(round(genome_label$chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.15, just = c("left", "bottom"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor))
      grid.text(label = paste(round(genome_label$chromendlabel, 1), scale, sep = " "), x = 1, y = 0.15, just = c("right", "bottom"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor))

    } else if (side == "left"){

      grid.segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, gp = gpar(col = genome_label$linecolor, lwd = genome_label$lwd))
      grid.text(label = genome_label$chromlabel, x = 0.85, y = 0.5, gp = gpar(fontface = "bold", fontsize = genome_label$fontsize, col = genome_label$fontcolor), just = c("center", "bottom"), rot = 90)
      grid.text(label = paste(round(genome_label$chromstartlabel, 1), scale, sep = " "), x = 0.85, y = 0, just = c("left", "bottom"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor), rot = 90)
      grid.text(label = paste(round(genome_label$chromendlabel, 1), scale, sep = " "), x = 0.85, y = 1, just = c("right", "bottom"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor), rot = 90)


    } else if (side == "right"){

      grid.segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, gp = gpar(col = genome_label$linecolor, lwd = genome_label$lwd))
      grid.text(label = genome_label$chromlabel, x = 0.15, y = 0.5, gp = gpar(fontface = "bold", fontsize = genome_label$fontsize, col = genome_label$fontcolor), just = c("center", "top"), rot = 90)
      grid.text(label = paste(round(genome_label$chromstartlabel, 1), scale, sep = " "), x = 0.15, y = 0, just = c("left", "top"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor), rot = 90)
      grid.text(label = paste(round(genome_label$chromendlabel, 1), scale, sep = " "), x = 0.15, y = 1, just = c("right", "top"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor), rot = 90)

    }

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  genome_label <- structure(list(scale = scale, fontsize = fontsize, fontcolor = fontcolor, linecolor = linecolor, lwd = lwd), class = "genome_label")

  # ======================================================================================================================================================================================
  # SET UP PAGE/SCALE
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

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

  tg <- textGrob(label = scale, x = 0.5, y = 0.5, default.units = "npc", gp = gpar(fontsize = fontsize))

  # ======================================================================================================================================================================================
  # INPUT PLOT
  # ======================================================================================================================================================================================

  if (!is.null(plot)){

    chrom <- paste0("chr", plot$chrom)
    chromstartlabel <- plot$chromstart/fact
    chromendlabel <- plot$chromend/fact

    genome_label[["chromlabel"]] = chrom
    genome_label[["chromstartlabel"]] = chromstartlabel
    genome_label[["chromendlabel"]] = chromendlabel

    adjusted_coords <- adjust_coords(plot = plot, page_units = page_units, page_height = page_height)

    message(paste("Plot to label detected."))

    side <- readline(prompt = "Enter side to add genome label. Options are 'bottom', 'top', 'left', or 'right': ")

    genome_label <- label_coords(plot = plot, xCoord = adjusted_coords[[1]], yCoord = adjusted_coords[[2]], side = side, genome_label = genome_label)

    genome_label <- label_dims(plot = plot, side = side, text = tg, genome_label = genome_label)

    vp <- viewport(width = convertWidth(unit(genome_label$width, units = genome_label$units), unitTo = page_units),
                   height = convertHeight(unit(genome_label$height, units = genome_label$units), unitTo = page_units),
                   x = convertX(unit(genome_label$x, units = genome_label$units), unitTo = page_units),
                   y = convertY(unit(genome_label$y, units = genome_label$units), unitTo = page_units), just = genome_label$just)
    pushViewport(vp)

    plot_label(side = side, genome_label = genome_label)

  # ======================================================================================================================================================================================
  # CUSTOMIZED PLOT
  # ======================================================================================================================================================================================

  } else {

    chromstartlabel = chromstart/fact
    chromendlabel = chromend/fact

    label_width <- convertWidth(unit(width, units = units), unitTo = page_units, valueOnly = TRUE)
    label_height <- convertHeight(heightDetails(tg), unitTo = page_units, valueOnly = TRUE)
    new_x <- convertX(unit(x, unit = units), unitTo = page_units, valueOnly = TRUE)
    new_y <- convertY(unit(y, unit = units), unitTo = page_units, valueOnly = TRUE)

    genome_label[["chromlabel"]] = chrom
    genome_label[["chromstartlabel"]] = chromstartlabel
    genome_label[["chromendlabel"]] = chromendlabel
    genome_label[["x"]] = x
    genome_label[["y"]] = y
    genome_label[["units"]] = units
    genome_label[["just"]] = just
    genome_label[["width"]] = width
    genome_label[["height"]] = convertHeight(heightDetails(tg), unitTo = units, valueOnly = TRUE)

    vp <- viewport(width = unit(label_width, page_units), height = unit(label_height, page_units),
                   x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), just = just, angle = rotation)
    pushViewport(vp)


    if (lineAbove == T){

      grid.segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = linecolor, lwd = lwd))
      grid.text(label = chrom, x = 0.5, y = 0.85, gp = gpar(fontface = "bold", fontsize = fontsize, col = fontcolor), just = c("center", "top"))
      grid.text(label = paste(round(chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.85, just = c("left", "top"), gp = gpar(fontsize = fontsize, col = linecolor))
      grid.text(label = paste(round(chromendlabel, 1), scale, sep = " "), x = 1, y = 0.85, just = c("right","top"), gp = gpar(fontsize = fontsize, col = linecolor))

    } else if (lineAbove == F){

      grid.segments(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(col = linecolor, lwd = lwd))
      grid.text(label = chrom, x = 0.5, y = 0.15, gp = gpar(fontface = "bold", fontsize = fontsize, col = fontcolor), just = c("center", "bottom"))
      grid.text(label = paste(round(chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.15, just = c("left", "bottom"), gp = gpar(fontsize = fontsize, col = linecolor))
      grid.text(label = paste(round(chromendlabel, 1), scale, sep = " "), x = 1, y = 0.15, just = c("right", "bottom"), gp = gpar(fontsize = fontsize, col = linecolor))

    }


  }

  ## Go back up a viewport
  upViewport()
}
