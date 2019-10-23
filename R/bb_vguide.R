bb_vguide <- function(y, units="inches", col="grey55", envir=bbEnv){
  g <- segmentsGrob(x0=y, x1=y, y0=unit(0, units = "npc"), y1=unit(1, units = "npc"), default.units = units, gp = gpar(col = col))
  appendGuide(g)
}
