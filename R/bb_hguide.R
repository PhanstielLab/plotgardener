bb_hguide <- function(x, units="inches", col="grey55", envir=bbEnv){
  g <- segmentsGrob(x0=unit(0, units = "npc"), x1=unit(1, units = "npc"), y0=envir$page_height-x, y1=envir$page_height-x, default.units = units,
                    gp = gpar(col = col))
  appendGuide(g)
}
