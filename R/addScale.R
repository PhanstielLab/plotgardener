addScale <- function(colVec, minLab, maxLab, x=0.93, y=0.72, width=0.03, height=0.3, just=c("center", "center"),
                     units="npc", fontsize=10, fontcolor="grey75", fontface = "plain", orientation = "vertical", border = T){

  plotArea <- viewport(x=x, y=y, width = width, height = height, just = just, default.units = units)
  scaleArea <- viewport(x=0.5, y=0.5, width = 1, height = 1)

  pushViewport(plotArea)
  pushViewport(scaleArea)

  if (orientation == "vertical"){
    lowLab <- textGrob(label = minLab, x = 0.5, y = 0, just = c("center", "top"), default.units = "npc", gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    highLab <- textGrob(label = maxLab, x = 0.5, y = 1, just = c("center", "bottom"), default.units = "npc", gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    lH <- convertHeight(x = heightDetails(lowLab), unitTo = "npc", valueOnly = T)
    hH <- convertHeight(x = heightDetails(highLab), unitTo = "npc", valueOnly = T)

    grid.draw(lowLab)
    grid.raster(rev(colVec), x = 0.5, y = 0.5, width = 1, height = unit(1-lH, "npc")) # for some reason you only need lH or hH but not both
    grid.draw(highLab)
    if (border == T) {
      grid.rect(x = 0.5, y = 0.5, width = 1, height = unit(1-lH, "npc"), gp = gpar(col = "black", fill = NA))
    }

  } else if (orientation == "horizontal"){
    lowLab <- textGrob(label = minLab, x = 0, y = 0.5, just = c("left", "center"), default.units = "npc", gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    highLab <- textGrob(label = maxLab, x = 1, y = 0.5, just = c("right", "center"), default.units = "npc", gp = gpar(fontsize = fontsize, col = fontcolor, fontface = fontface))
    lW <- convertWidth(x = widthDetails(lowLab), unitTo = "npc")
    hW <- convertWidth(x = widthDetails(highLab), unitTo = "npc")

    grid.draw(lowLab)
    grid.raster(matrix(data = colVec, nrow = 1, ncol = length(colVec)), x = unit(0, "npc")+lW, y = 0.5, width = (unit(1, "npc")-hW-lW*(1/2)), height = 1, just = c("left", "center"))
    grid.draw(highLab)
    if (border == T) {
      grid.rect(x = unit(0, "npc")+lW, y = 0.5, width = (unit(1, "npc")-hW-lW*(1/2)), height = 1, just = c("left", "center"), gp = gpar(col = "black", fill = NA))
    }
  }

  upViewport()
  upViewport()
}
