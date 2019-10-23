bb_makepage <- function(width = 8.5, height = 11, units = "inches", showOutline = TRUE, showRuler = TRUE, xgrid = 0.5, ygrid = 0.5){

  grid.newpage()
  bbEnv$guides <- list()

  page_vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                      width = unit(width, units = units), height = unit(height, units = units),
                      xscale = c(0, width),
                      name = "bb_page")

  pushViewport(page_vp)

  if (showOutline == TRUE){
    border <- rectGrob()
    appendGuide(border)

  }

  ## Draw ruler segments
  if (showRuler == TRUE){
    div <- 1/16
    x0 <- 0
    x1 <- -1/32
    y0 <- height
    y1 <- height + 1/32
    xsegs <- segmentsGrob(x0 = seq(0, width, div), y0 = y0,
                          x1 = seq(0, width, div), y1 = y1, default.units = units, gp = gpar(col = "black"))
    ysegs <- segmentsGrob(x0 = x0, y0 = seq(0, height, div),
                          x1 = x1, y1 = seq(0, height, div), default.units = units, gp = gpar(col = "black"))
    # grid.draw(xsegs)
    # grid.draw(ysegs)
    appendGuide(xsegs)
    appendGuide(ysegs)
    for (i in 1:4){
      div <- div*2
      x0 <- x1
      x1 <- x1 - 1/32
      y0 <- y1
      y1 <- y1 + 1/32
      v <- segmentsGrob(x0 = xsegs$x0[xsegs$x0 %in% seq(0, width, div)], y0 = y0,
                        x1 = xsegs$x0[xsegs$x0 %in% seq(0, width, div)], y1 = y1, default.units = units, gp = gpar(col = "black"))
      appendGuide(v)
      h <- segmentsGrob(x0 = x0, y0 = ysegs$y1[ysegs$y1 %in% seq(0, height, div)],
                        x1 = x1, y1 = ysegs$y1[ysegs$y1 %in% seq(0, height, div)], default.units = units, gp = gpar(col = "black"))
      appendGuide(h)
    }
    hLabel <- textGrob(label = seq(0, width, div), x = seq(0, width, div), y = y0, vjust = -0.5, default.units = units)
    appendGuide(hLabel)
    vLabel <- textGrob(label = seq(0, height, div), x = x0, y = rev(seq(0, height, div)), hjust = 1.5, default.units = units)
    appendGuide(vLabel)

    ## Unit annotation
    unitLabel <- textGrob(label = units, x = 0, y = height, hjust = 1.75, vjust = -1.5, default.units = units, just = c("right", "bottom"))
    appendGuide(unitLabel)
  }


  tryCatch(expr = {
    ## Draw x and y gridlines
    xGrid <- segmentsGrob(x0 = seq(0, width, xgrid), y0 = 0,
                          x1 = seq(0, width, xgrid), y1 = height,
                          default.units = units, gp = gpar(col = "grey50", lty = 2, lwd = 0.5))
    appendGuide(xGrid)

    yGrid <- segmentsGrob(x0 = 0, y0 = seq(0, height, ygrid),
                          x1 = width, y1 = seq(0, height, ygrid),
                          default.units = units, gp = gpar(col = "grey50", lty = 2, lwd = 0.5))
    appendGuide(yGrid)

  }, error = function(e) return())




  assign("page_height", height, envir = bbEnv)
  assign("page_units", units, envir = bbEnv)

}
