#' Adds genome coordinates to the axis of a plot
#'
#' @param plot plot to annotate (saved as output from BentoBox plotting function)
#' @param chrom if not inputting plot, chromosome of region to plot
#' @param chromstart if not inputting plot, chromstart of region to plot
#' @param chromend if not inputting plot, chromend of region to plot
#' @param width if not inputting plot, width of label
#' @param x if not inputting plot, x-coordinate of label
#' @param y if not inputting plot, y-coordinate of label
#' @param just justification
#' @param units units of inputted x and y
#' @param scale scale of the plot; options are "bp", "Kb", or "Mb"
#' @param side side of the plot to annotate, options are "bottom", "top", "right", or "left"
#' @param fontsize fontsize for labels (in points)
#' @param fontcolor color for all text/lines
#' @param lineAbove TRUE/FALSE indicating if line should be above or below text
#' @param linecolor linecolor
#' @param lwd linewidth


#' @export

bb_labelGenome <- function(plot = NULL, chrom, chromstart, chromend, width, x, y, just = c("left", "top"), rotation = 0,
                           units = "inches", scale = "bp", fontsize = 10, fontcolor = "black", lineAbove = T, linecolor = "black", lwd = 1){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_labelgenome
  # errorcheck_bb_labelgenome <- function()



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

      line <- grid.segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = genome_label$linecolor, lwd = genome_label$lwd))
      chromLab <- grid.text(label = genome_label$chromlabel, x = 0.5, y = 0.85, gp = gpar(fontface = "bold", fontsize = genome_label$fontsize, col = genome_label$fontcolor), just = c("center", "top"))
      startLab <- grid.text(label = paste(round(genome_label$chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.85, just = c("left", "top"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor))
      endLab <- grid.text(label = paste(round(genome_label$chromendlabel, 1), scale, sep = " "), x = 1, y = 0.85, just = c("right","top"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor))

    } else if (side == "top"){

      line <- grid.segments(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(col = genome_label$linecolor, lwd = genome_label$lwd))
      chromLab <- grid.text(label = genome_label$chromlabel, x = 0.5, y = 0.15, gp = gpar(fontface = "bold", fontsize = genome_label$fontsize, col = genome_label$fontcolor), just = c("center", "bottom"))
      startLab <- grid.text(label = paste(round(genome_label$chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.15, just = c("left", "bottom"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor))
      endLab <- grid.text(label = paste(round(genome_label$chromendlabel, 1), scale, sep = " "), x = 1, y = 0.15, just = c("right", "bottom"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor))

    } else if (side == "left"){

      line <- grid.segments(x0 = 1, x1 = 1, y0 = 0, y1 = 1, gp = gpar(col = genome_label$linecolor, lwd = genome_label$lwd))
      chromLab <- grid.text(label = genome_label$chromlabel, x = 0.85, y = 0.5, gp = gpar(fontface = "bold", fontsize = genome_label$fontsize, col = genome_label$fontcolor), just = c("center", "bottom"), rot = 90)
      startLab <- grid.text(label = paste(round(genome_label$chromstartlabel, 1), scale, sep = " "), x = 0.85, y = 0, just = c("left", "bottom"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor), rot = 90)
      endLab <- grid.text(label = paste(round(genome_label$chromendlabel, 1), scale, sep = " "), x = 0.85, y = 1, just = c("right", "bottom"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor), rot = 90)

    } else if (side == "right"){

      line <- grid.segments(x0 = 0, x1 = 0, y0 = 0, y1 = 1, gp = gpar(col = genome_label$linecolor, lwd = genome_label$lwd))
      chromLab <- grid.text(label = genome_label$chromlabel, x = 0.15, y = 0.5, gp = gpar(fontface = "bold", fontsize = genome_label$fontsize, col = genome_label$fontcolor), just = c("center", "top"), rot = 90)
      startLab <- grid.text(label = paste(round(genome_label$chromstartlabel, 1), scale, sep = " "), x = 0.15, y = 0, just = c("left", "top"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor), rot = 90)
      endLab <- grid.text(label = paste(round(genome_label$chromendlabel, 1), scale, sep = " "), x = 0.15, y = 1, just = c("right", "top"), gp = gpar(fontsize = genome_label$fontsize, col = genome_label$fontcolor), rot = 90)

    }

    assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_genome_label <- structure(list(scale = scale, grobs = NULL, viewport = NULL,
                                    gpar = list(fontsize = fontsize, fontcolor = fontcolor, linecolor = linecolor, lwd = lwd)), class = "genome_label")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage()

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
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("label_grobs", gTree(name = "label_grobs"), envir = bbEnv)

  # ======================================================================================================================================================================================
  # INPUT PLOT
  # ======================================================================================================================================================================================

  if (!is.null(plot)){

    chrom <- paste0("chr", plot$chrom)
    chromstartlabel <- plot$chromstart/fact
    chromendlabel <- plot$chromend/fact

    bb_genome_label[["chromlabel"]] = chrom
    bb_genome_label[["chromstartlabel"]] = chromstartlabel
    bb_genome_label[["chromendlabel"]] = chromendlabel

    adjusted_coords <- adjust_coords(plot = plot, page_units = page_units, page_height = page_height)

    message(paste(class(plot), "plot to label detected."))

    side <- readline(prompt = "Enter side to add genome label. Options are 'bottom', 'top', 'left', or 'right': ")

    bb_genome_label <- label_coords(plot = plot, xCoord = adjusted_coords[[1]], yCoord = adjusted_coords[[2]], side = side, genome_label = bb_genome_label)

    bb_genome_label <- label_dims(plot = plot, side = side, text = tg, genome_label = bb_genome_label)

    ## Name viewport
    current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
    vp_name <- paste0("bb_genomeLabel", length(grep(pattern = "bb_genomeLabel", x = current_viewports)) + 1)

    vp <- viewport(width = convertWidth(unit(bb_genome_label$width, units = bb_genome_label$units), unitTo = page_units),
                   height = convertHeight(unit(bb_genome_label$height, units = bb_genome_label$units), unitTo = page_units),
                   x = convertX(unit(bb_genome_label$x, units = bb_genome_label$units), unitTo = page_units),
                   y = convertY(unit(bb_genome_label$y, units = bb_genome_label$units), unitTo = page_units), just = bb_genome_label$just, name = vp_name)
    bb_genome_label$viewport <- vp
    pushViewport(vp)

    plot_label(side = side, genome_label = bb_genome_label)

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

    bb_genome_label[["chromlabel"]] = chrom
    bb_genome_label[["chromstartlabel"]] = chromstartlabel
    bb_genome_label[["chromendlabel"]] = chromendlabel
    bb_genome_label[["x"]] = x
    bb_genome_label[["y"]] = y
    bb_genome_label[["units"]] = units
    bb_genome_label[["just"]] = just
    bb_genome_label[["width"]] = width
    bb_genome_label[["height"]] = convertHeight(heightDetails(tg), unitTo = units, valueOnly = TRUE)

    ## Name viewport
    current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
    vp_name <- paste0("bb_genomeLabel", length(grep(pattern = "bb_genomeLabel", x = current_viewports)) + 1)

    vp <- viewport(width = unit(label_width, page_units), height = unit(label_height, page_units),
                   x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), just = just, angle = rotation, name = vp_name)
    bb_genome_label$viewport <- vp
    pushViewport(vp)

    if (lineAbove == T){

      line <- grid.segments(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = linecolor, lwd = lwd))
      chromLab <- grid.text(label = chrom, x = 0.5, y = 0.85, gp = gpar(fontface = "bold", fontsize = fontsize, col = fontcolor), just = c("center", "top"))
      startLab <- grid.text(label = paste(round(chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.85, just = c("left", "top"), gp = gpar(fontsize = fontsize, col = linecolor))
      endLab <- grid.text(label = paste(round(chromendlabel, 1), scale, sep = " "), x = 1, y = 0.85, just = c("right","top"), gp = gpar(fontsize = fontsize, col = linecolor))

    } else if (lineAbove == F){

      line <- grid.segments(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(col = linecolor, lwd = lwd))
      chromLab <- grid.text(label = chrom, x = 0.5, y = 0.15, gp = gpar(fontface = "bold", fontsize = fontsize, col = fontcolor), just = c("center", "bottom"))
      startLab <- grid.text(label = paste(round(chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.15, just = c("left", "bottom"), gp = gpar(fontsize = fontsize, col = linecolor))
      endLab <- grid.text(label = paste(round(chromendlabel, 1), scale, sep = " "), x = 1, y = 0.15, just = c("right", "bottom"), gp = gpar(fontsize = fontsize, col = linecolor))

    }

    assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)
  }

  ## Go back up a viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ASSIGN GROBS TO SCALE OBJECT
  # ======================================================================================================================================================================================

  ## Add grobs to scale object
  bb_genome_label$grobs <- get("label_grobs", envir = bbEnv)$children

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_genome_label)
}
