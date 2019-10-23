## Function to add genome label ####
bb_labelGenome <- function(object = mth, x, y, just=c("left", "top"), units="inches",
                           unitLabel = "Mb", lwd = 4, color = "grey75", fontsize = 10, fontface = "bold"){

  ## Convert and reassign coordinates ####
  # convertedCoordinates <- bb_convertCoordinates(x, y, width = object$width, height = 0.2, units = units)
  # x <- convertedCoordinates$x
  # y <- convertedCoordinates$y
  # width <- convertedCoordinates$width
  # height <- convertedCoordinates$height

  ## Assign units to scale ####
  if (unitLabel == "Mb"){
    sizeFactor <- 10^6
  } else if (unitLabel == "Kb"){
    sizeFactor <- 10^3
  } else if (unitLabel == "bp"){
    sizeFactor <- 1
  }

  ## Define string height for viewpoint ####
  tg <- textGrob(label = unitLabel, x = 0.5, y = 0.5, default.units = "npc", gp = gpar(fontsize = fontsize))

  ## Define and push viewport ####
  vp <- viewport(x = unit(x, units), y = unit(y, units), just = just,
                 width = unit(object$width, object$units), height = heightDetails(tg)*fontsize*0.1)
  pushViewport(vp)

  # grid.rect(gp = gpar(col = "black"))

  ## Axis line ####
  grid.segments(x0 = 0, y0 = 1, x1 = 1, y1 = 1, default.units = "npc", gp = gpar(lwd = lwd, col = color, lineend = 1))

  ## Axis labels ####
  grid.text(label = paste(round(object$chromstart/sizeFactor, 1), unitLabel),
            x = 0, y = 0.85, just = c("left", "top"), gp = gpar(fontsize = fontsize, fontface = fontface, col = color))
  grid.text(label = paste(round(object$chromend/sizeFactor, 1), unitLabel),
            x = 1, y = 0.85, just = c("right", "top"), gp = gpar(fontsize = fontsize, fontface = fontface, col = color))
  grid.text(label = object$chrom, x = 0.5, y = 0.85, just = c("center", "top"), gp = gpar(fontsize = fontsize, fontface = fontface, col = color))

  ## Exit viewport ####
  upViewport()
}
