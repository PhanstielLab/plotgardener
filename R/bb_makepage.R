#' Makes a page for a BentoBox layout
#'
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param width A numeric specifying width
#' @param height A numeric ospecifying height
#' @param default.units A string indicating the units of width and height
#' @param showOutline TRUE/FALSE indicating whether to draw a black outline around the entire page
#' @param showRuler TRUE/FALSE indicating whether to show guiding ruler along top and side of page
#' @param xgrid vertical gridlines
#' @param ygrid horizontal gridlines
#' @export

bb_makePage <- function(params = NULL, width = 8.5, height = 11, default.units = "inches", showOutline = TRUE, showRuler = TRUE, xgrid = 0.5, ygrid = 0.5){

  # ======================================================================================================================================================================================
  # MAKE NEW PAGE
  # ======================================================================================================================================================================================

  grid.newpage()

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(width)) width <- NULL
  if(missing(height)) height <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(showOutline)) showOutline <- NULL
  if(missing(showRuler)) showRuler <- NULL
  if(missing(xgrid)) xgrid <- NULL
  if(missing(ygrid)) ygrid <- NULL

  ## Compile all parameters into an internal object
  bb_page <- structure(list(width = width, height = height, default.units = default.units, showOutline = showOutline,
                            showRuler = showRuler, xgrid = xgrid, ygrid = ygrid), class = "bb_page")
  bb_page <- parseParams(bb_params = params, object_params = bb_page)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_page$width)) bb_page$width <- 8.5
  if(is.null(bb_page$height)) bb_page$height <- 11
  if(is.null(bb_page$default.units)) bb_page$default.units <- "inches"
  if(is.null(bb_page$showOutline)) bb_page$showOutline <- TRUE
  if(is.null(bb_page$showRuler)) bb_page$showRuler <- TRUE
  if(is.null(bb_page$xgrid)) bb_page$xgrid <- 0.5
  if(is.null(bb_page$ygrid)) bb_page$ygrid <- 0.5

  # ======================================================================================================================================================================================
  # VIEWPORT
  # =====================================================================================================================================================================================

  page_vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                      width = unit(bb_page$width, units = bb_page$default.units), height = unit(bb_page$height, units = bb_page$default.units),
                      xscale = c(0, bb_page$width), yscale = rev(c(0, bb_page$height)),
                      name = "bb_page")

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("guide_grobs", gTree(name = "guide_grobs", vp = page_vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # SHOW OUTLINE
  # ======================================================================================================================================================================================

  if (bb_page$showOutline == TRUE){

    border <- rectGrob()
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = border), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # SHOW RULER
  # ======================================================================================================================================================================================

  if (bb_page$showRuler == TRUE){

    div <- 1/16
    x0 <- 0
    x1 <- -1/32
    y0 <- bb_page$height
    y1 <- bb_page$height + 1/32
    xsegs <- segmentsGrob(x0 = seq(0, bb_page$width, div), y0 = y0,
                          x1 = seq(0, bb_page$width, div), y1 = y1, default.units = bb_page$default.units, gp = gpar(col = "black"))
    ysegs <- segmentsGrob(x0 = x0, y0 = seq(0, bb_page$height, div),
                          x1 = x1, y1 = seq(0, bb_page$height, div), default.units = bb_page$default.units, gp = gpar(col = "black"))
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = xsegs), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = ysegs), envir = bbEnv)

    for (i in 1:4){
      div <- div*2
      x0 <- x1
      x1 <- x1 - 1/32
      y0 <- y1
      y1 <- y1 + 1/32
      v <- segmentsGrob(x0 = xsegs$x0[xsegs$x0 %in% seq(0, bb_page$width, div)], y0 = y0,
                        x1 = xsegs$x0[xsegs$x0 %in% seq(0, bb_page$width, div)], y1 = y1, default.units = bb_page$default.units, gp = gpar(col = "black"))
      h <- segmentsGrob(x0 = x0, y0 = ysegs$y1[ysegs$y1 %in% seq(0, bb_page$height, div)],
                        x1 = x1, y1 = ysegs$y1[ysegs$y1 %in% seq(0, bb_page$height, div)], default.units = bb_page$default.units, gp = gpar(col = "black"))

      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v), envir = bbEnv)
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h), envir = bbEnv)
    }

    hLabel <- textGrob(label = seq(0, bb_page$width, div), x = seq(0, bb_page$width, div), y = y0, vjust = -0.5, default.units = bb_page$default.units)
    vLabel <- textGrob(label = seq(0, bb_page$height, div), x = x0, y = seq(0, bb_page$height, div), hjust = 1.5, default.units = "native")

    ## Unit annotation
    unitLabel <- textGrob(label = bb_page$default.units, x = 0, y = bb_page$height, hjust = 1.75, vjust = -1.5, default.units = bb_page$default.units, just = c("right", "bottom"))

    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = hLabel), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = vLabel), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = unitLabel), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # X AND Y GRIDLINES
  # ======================================================================================================================================================================================

  ## Draw x and y gridlines
  tryCatch(expr = {

    xGrid <- segmentsGrob(x0 = seq(0, bb_page$width, bb_page$xgrid), y0 = 0,
                          x1 = seq(0, bb_page$width, bb_page$xgrid), y1 = bb_page$height,
                          default.units = bb_page$default.units, gp = gpar(col = "grey50", lty = 2, lwd = 0.5))

    yGrid <- segmentsGrob(x0 = 0, y0 = seq(0, bb_page$height, bb_page$ygrid),
                          x1 = bb_page$width, y1 = seq(0, bb_page$height, bb_page$ygrid),
                          default.units = bb_page$default.units, gp = gpar(col = "grey50", lty = 2, lwd = 0.5))

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

  assign("page_height", bb_page$height, envir = bbEnv)
  assign("page_units", bb_page$default.units, envir = bbEnv)

}
