addApaScale <- function(apaObj, x=0.93, y=0.72, width=0.03, height=0.3, just=c("center", "center"), units="npc", fontsize=10, fontcolor="grey75"){

  plotArea <- viewport(x=x, y=y, width = width, height = height, just = just, default.units = units)
  scaleArea <- viewport(x=0.5, y=0.5, width = 1, height = 1)

  pushViewport(plotArea)
  pushViewport(scaleArea)

  grid.raster(rev(apaObj$col_fun(1000)), x = 0.5, y = 0.5, width = 1, height = 1-(2.5*fontsize/100), default.units = "npc")

  if(is.null(apaObj$zrange)){
    grid.text(label = round(max(apaObj$dat)), x = unit(x = 0.5, units = "npc"), y = unit(x = 1, units = "npc"),
              gp = gpar(fontsize = fontsize, fontface = "bold", col = fontcolor), just = c("center", "top"))
    grid.text(label = round(min(apaObj$dat)), x = unit(x = 0.5, units = "npc"), y = unit(x = 0, units = "npc"),
              gp = gpar(fontsize = fontsize, fontface = "bold", col = fontcolor), just = c("center", "bottom"))
  } else {
    grid.text(label = round(apaObj$zrange[2]), x = unit(x = 0.5, units = "npc"), y = unit(x = 1, units = "npc"),
              gp = gpar(fontsize = fontsize, fontface = "bold", col = fontcolor), just = c("center", "top"))
    grid.text(label = round(apaObj$zrange[1]), x = unit(x = 0.5, units = "npc"), y = unit(x = 0, units = "npc"),
              gp = gpar(fontsize = fontsize, fontface = "bold", col = fontcolor), just = c("center", "bottom"))
  }

  upViewport()
  upViewport()
}
