arrowAnnotation <- function(hic = sth, lineDist = 0.25, shift = 10000, loop = loop1){
  pushViewport(viewport(x = hic$x, y = height-hic$y, width = hic$width, height = hic$height, default.units = "inches", just = c("left", "top"),
                        xscale = c(hic$chromstart, hic$chromend), yscale = c(hic$chromstart, hic$chromend)))
  grid.lines(x = unit(c(loop$x1-span1*lineDist, loop$x1-shift), "native"),
             y = unit(c(loop$y1+span1*lineDist, loop$y1+shift), "native"),
             gp = gpar(fill="black"),
             arrow = arrow(length = unit(0.1, "inches"),
                           ends="last", type="closed"))
  grid.lines(x = unit(c(loop$y1+span1*lineDist, loop$y1+shift), "native"),
             y = unit(c(loop$x1-span1*lineDist, loop$x1-shift), "native"),
             gp = gpar(fill="black"),
             arrow = arrow(length = unit(0.1, "inches"),
                           ends="last", type="closed"))
  upViewport()
}
