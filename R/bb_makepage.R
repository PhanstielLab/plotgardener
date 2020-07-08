#' Makes a page for a BentoBox layout
#'
#' @param width A numeric specifying width
#' @param height A numeric ospecifying height
#' @param default.units A string indicating the units of width and height
#' @param showOutline TRUE/FALSE indicating whether to draw a black outline around the entire page
#' @param showRuler TRUE/FALSE indicating whether to show guiding ruler along top and side of page
#' @param xgrid vertical gridlines
#' @param ygrid horizontal gridlines
#' @export

bb_makePage <- function(width = 8.5, height = 11, default.units = "inches", showOutline = TRUE, showRuler = TRUE, xgrid = 0.5, ygrid = 0.5){

  # ======================================================================================================================================================================================
  # MAKE NEW PAGE
  # ======================================================================================================================================================================================

  grid.newpage()

  # ======================================================================================================================================================================================
  # VIEWPORT
  # =====================================================================================================================================================================================

  page_vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                      width = unit(width, units = default.units), height = unit(height, units = default.units),
                      xscale = c(0, width), yscale = rev(c(0, height)),
                      name = "bb_page")

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("guide_grobs", gTree(name = "guide_grobs", vp = page_vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # SHOW OUTLINE
  # ======================================================================================================================================================================================

  if (showOutline == TRUE){

    border <- rectGrob()
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = border), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # SHOW RULER
  # ======================================================================================================================================================================================

  if (showRuler == TRUE){

    div <- 1/16
    x0 <- 0
    x1 <- -1/32
    y0 <- height
    y1 <- height + 1/32
    xsegs <- segmentsGrob(x0 = seq(0, width, div), y0 = y0,
                          x1 = seq(0, width, div), y1 = y1, default.units = default.units, gp = gpar(col = "black"))
    ysegs <- segmentsGrob(x0 = x0, y0 = seq(0, height, div),
                          x1 = x1, y1 = seq(0, height, div), default.units = default.units, gp = gpar(col = "black"))
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = xsegs), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = ysegs), envir = bbEnv)

    for (i in 1:4){
      div <- div*2
      x0 <- x1
      x1 <- x1 - 1/32
      y0 <- y1
      y1 <- y1 + 1/32
      v <- segmentsGrob(x0 = xsegs$x0[xsegs$x0 %in% seq(0, width, div)], y0 = y0,
                        x1 = xsegs$x0[xsegs$x0 %in% seq(0, width, div)], y1 = y1, default.units = default.units, gp = gpar(col = "black"))
      h <- segmentsGrob(x0 = x0, y0 = ysegs$y1[ysegs$y1 %in% seq(0, height, div)],
                        x1 = x1, y1 = ysegs$y1[ysegs$y1 %in% seq(0, height, div)], default.units = default.units, gp = gpar(col = "black"))

      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v), envir = bbEnv)
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h), envir = bbEnv)
    }

    hLabel <- textGrob(label = seq(0, width, div), x = seq(0, width, div), y = y0, vjust = -0.5, default.units = default.units)
    vLabel <- textGrob(label = seq(0, height, div), x = x0, y = seq(0, height, div), hjust = 1.5, default.units = "native")

    ## Unit annotation
    unitLabel <- textGrob(label = default.units, x = 0, y = height, hjust = 1.75, vjust = -1.5, default.units = default.units, just = c("right", "bottom"))

    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = hLabel), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = vLabel), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = unitLabel), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # X AND Y GRIDLINES
  # ======================================================================================================================================================================================

  ## Draw x and y gridlines
  tryCatch(expr = {

    xGrid <- segmentsGrob(x0 = seq(0, width, xgrid), y0 = 0,
                          x1 = seq(0, width, xgrid), y1 = height,
                          default.units = default.units, gp = gpar(col = "grey50", lty = 2, lwd = 0.5))

    yGrid <- segmentsGrob(x0 = 0, y0 = seq(0, height, ygrid),
                          x1 = width, y1 = seq(0, height, ygrid),
                          default.units = default.units, gp = gpar(col = "grey50", lty = 2, lwd = 0.5))

    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = xGrid), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = yGrid), envir = bbEnv)

  }, error = function(e) return())

  # ======================================================================================================================================================================================
  # DRAW GROBS
  # ======================================================================================================================================================================================

  grid.draw(get("guide_grobs", envir = bbEnv))
  downViewport("bb_page")

  # ======================================================================================================================================================================================
  # ASSIGN PAGE PARAMETERS
  # ======================================================================================================================================================================================

  assign("page_height", height, envir = bbEnv)
  assign("page_units", default.units, envir = bbEnv)

}
