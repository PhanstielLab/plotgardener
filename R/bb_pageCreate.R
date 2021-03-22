#' Create a page for a BentoBox layout
#'
#' @param width A numeric or unit object specifying page width. Default value is \code{width = 8}.
#' @param height A numeric or unit object specifying page height. Default value is \code{height = 11}.
#' @param default.units A string indicating the default units to use if \code{width} or \code{height} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param xgrid A numeric indicating the increment by which to place vertical gridlines. Default value is \code{xgrid = 0.5}.
#' @param ygrid A numeric indicating the increment by which to place horizontal gridlines. Default value is \code{ygrid = 0.5}.
#' @param showGuides A logical value indicating whether to draw a black border around the entire page and guiding rulers along the top and left side of the page. Default value is \code{showOutline = TRUE}.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#'
#' @examples
#' ## Create a 6-inch wide, 4.5-inch high BentoBox page
#' bb_pageCreate(width = 6, height = 4.5, default.units = "inches")
#'
#' ## Create a 14-cm wide, 10-cm high Bentobox page
#' bb_pageCreate(width = 14, height = 10, default.units = "cm")
#'
#' @details \code{width} and \code{height} must be specified in the same units.
#'
#' @export
bb_pageCreate <- function(width = 8.5, height = 11, default.units = "inches", xgrid = 0.5, ygrid = 0.5, showGuides = TRUE, params = NULL){

  # ======================================================================================================================================================================================
  # FUNCTION
  # ======================================================================================================================================================================================

  errorcheck_bb_pageCreate <- function(object){

    ## Width and height can't be negative
    if (as.numeric(object$width) <= 0){
      stop("`width` cannot be 0 or negative.", call. = FALSE)
    }

    if (as.numeric(object$height) <= 0){
      stop("`height` cannot be 0 or negative.", call. = FALSE)
    }


    ## Width and height need to be of the same units
    widthUnit <- gsub("[0-9]|[.]", "", object$width)
    heightUnit <- gsub("[0-9]|[.]", "", object$height)
    if (widthUnit != heightUnit){

      stop(paste("`width` and `height` must be in the same units. `width` detected as", paste0("`",widthUnit,"`"),
                 "and `height` detected as", paste0(paste0("`",heightUnit,"`"), ".")), call. = FALSE)
    }

    ## Page units must be a valid option
    page_units <- gsub("[0-9]|[.]", "", bb_page$width)
    if (!page_units %in% validUnits){
      stop(paste(c("Invalid page units. Options for page units are:", validUnits), collapse = " "), call. = FALSE)
    }



  }

  # ======================================================================================================================================================================================
  # MAKE NEW PAGE
  # ======================================================================================================================================================================================

  grid.newpage()
  assign("bb_vpTree", list(), envir = bbEnv)

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(width)) width <- NULL
  if(missing(height)) height <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(showGuides)) showGuides <- NULL
  if(missing(xgrid)) xgrid <- NULL
  if(missing(ygrid)) ygrid <- NULL

  ## Compile all parameters into an internal object
  bb_page <- structure(list(width = width, height = height, default.units = default.units, showGuides = showGuides,
                            xgrid = xgrid, ygrid = ygrid), class = "bb_page")
  bb_page <- parseParams(bb_params = params, object_params = bb_page)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_page$width)) bb_page$width <- 8.5
  if(is.null(bb_page$height)) bb_page$height <- 11
  if(is.null(bb_page$default.units)) bb_page$default.units <- "inches"
  if(is.null(bb_page$showGuides)) bb_page$showGuides <- TRUE
  if(is.null(bb_page$xgrid)) bb_page$xgrid <- 0.5
  if(is.null(bb_page$ygrid)) bb_page$ygrid <- 0.5

  # ======================================================================================================================================================================================
  # DEFAULT UNITS
  # =====================================================================================================================================================================================


  if (!"unit" %in% class(bb_page$width)){

    if (!is.numeric(bb_page$width)){

      stop("`width` is neither a unit object nor a numeric value. Cannot create BentoBox page.", call. = FALSE)

    }

    if (is.null(bb_page$default.units)){

      stop("`width` detected as numeric. `default.units` must be specified.", call. = FALSE)

    }

    if (!bb_page$default.units %in% validUnits){
      stop(paste(c("Invalid default units. Options for page units are:", validUnits), collapse = " "), call. = FALSE)
    }

    bb_page$width <- unit(bb_page$width, bb_page$default.units)

  }

  if (!"unit" %in% class(bb_page$height)){

    if (!is.numeric(bb_page$height)){

      stop("`height` is neither a unit object or a numeric value. Cannot create BentoBox page.", call. = FALSE)

    }

    if (is.null(bb_page$default.units)){

      stop("`height` detected as numeric. `default.units` must be specified.", call. = FALSE)

    }


    if (!bb_page$default.units %in% validUnits){
      stop(paste(c("Invalid default units. Options for page units are:", validUnits), collapse = " "), call. = FALSE)
    }

    bb_page$height <- unit(bb_page$height, bb_page$default.units)

  }

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # =====================================================================================================================================================================================

  errorcheck_bb_pageCreate(object = bb_page)

  # ======================================================================================================================================================================================
  # ASSIGN PAGE PARAMETERS
  # =====================================================================================================================================================================================
  page_height <- as.numeric(bb_page$height)
  assign("page_height", page_height, envir = bbEnv)

  page_units <- gsub("[0-9]|[.]", "", bb_page$width)
  assign("page_units", page_units, envir = bbEnv)

  page_width <- as.numeric(bb_page$width)

  # ======================================================================================================================================================================================
  # VIEWPORT
  # =====================================================================================================================================================================================

  page_vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                      width = bb_page$width, height = bb_page$height,
                      xscale = c(0, page_width), yscale = rev(c(0, page_height)),
                      name = "bb_page")

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("guide_grobs", gTree(name = "guide_grobs", vp = page_vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # SHOW OUTLINE/RULER/UNITS
  # ======================================================================================================================================================================================

  if (bb_page$showGuides == TRUE){

    border <- rectGrob()
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = border), envir = bbEnv)


    if (page_units == "inches"){

      div <- 1/16
      x0 <- 0
      x1 <- -1/32
      y0 <- page_height
      y1 <- page_height + 1/32
      xsegs <- segmentsGrob(x0 = seq(0, page_width, div), y0 = y0,
                            x1 = seq(0, page_width, div), y1 = y1, default.units = page_units, gp = gpar(col = "black"))
      ysegs <- segmentsGrob(x0 = x0, y0 = seq(0, page_height, div),
                            x1 = x1, y1 = seq(0, page_height, div), default.units = "native", gp = gpar(col = "black"))

      for (i in 1:4){
        div <- div*2
        x0 <- x1
        x1 <- x1 - 1/32
        y0 <- y1
        y1 <- y1 + 1/32
        v <- segmentsGrob(x0 = xsegs$x0[xsegs$x0 %in% seq(0, page_width, div)], y0 = y0,
                          x1 = xsegs$x0[xsegs$x0 %in% seq(0, page_width, div)], y1 = y1, default.units = page_units, gp = gpar(col = "black"))
        h <- segmentsGrob(x0 = x0, y0 = ysegs$y1[ysegs$y1 %in% seq(0, page_height, div)],
                          x1 = x1, y1 = ysegs$y1[ysegs$y1 %in% seq(0, page_height, div)], default.units = "native", gp = gpar(col = "black"))

        assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v), envir = bbEnv)
        assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h), envir = bbEnv)
      }

      hLabel <- textGrob(label = seq(0, page_width, div), x = seq(0, page_width, div), y = y0, vjust = -0.5, default.units = page_units)
      vLabel <- textGrob(label = seq(0, page_height, div), x = x0, y = seq(0, page_height, div), hjust = 1.5, default.units = "native")

    } else if (page_units == "cm"){

      tickH <- convertUnit(x = unit(1/32, "inches"), unitTo = "cm", valueOnly = T)
      div <- 1/10
      x0 <- 0
      x1 <- -3*tickH
      y0 <- page_height
      y1 <- page_height + 3*tickH
      xsegs <- segmentsGrob(x0 = seq(0, page_width, div), y0 = y0,
                            x1 = seq(0, page_width, div), y1 = y1, default.units = page_units, gp = gpar(col = "black"))
      ysegs <- segmentsGrob(x0 = x0, y0 = seq(0, page_height, div),
                            x1 = x1, y1 = seq(0, page_height, div), default.units = "native", gp = gpar(col = "black"))


      div2 <- 1/2
      v2 <- segmentsGrob(x0 = seq(0, page_width, div2), y0 = page_height + 3*tickH,
                        x1 = seq(0, page_width, div2), y1 = page_height + 4*tickH, default.units = page_units, gp = gpar(col = "black"))
      h2 <- segmentsGrob(x0 = -3*tickH, y0 = seq(0, page_height, div2),
                        x1 = -4*tickH, y1 = seq(0, page_height, div2), default.units = "native", gp = gpar(col = "black"))
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v2), envir = bbEnv)
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h2), envir = bbEnv)

      div3 <- 1
      v3 <- segmentsGrob(x0 = seq(0, page_width, div3), y0 = page_height + 4*tickH,
                         x1 = seq(0, page_width, div3), y1 = page_height + 5*tickH, default.units = page_units, gp = gpar(col = "black"))
      h3 <- segmentsGrob(x0 = -4*tickH, y0 = seq(0, page_height, div3),
                         x1 = -5*tickH, y1 = seq(0, page_height, div3), default.units = "native", gp = gpar(col = "black"))
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v3), envir = bbEnv)
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h3), envir = bbEnv)

      hLabel <- textGrob(label = seq(0, page_width, div3), x = seq(0, page_width, div3), y = page_height + 4*tickH, vjust = -0.5, default.units = page_units)
      vLabel <- textGrob(label = seq(0, page_height, div3), x = -4*tickH, y = seq(0, page_height, div3), hjust = 1.5, default.units = "native")


    } else if (page_units == "mm"){

      tickH <- convertUnit(x = unit(1/32, "inches"), unitTo = "mm", valueOnly = T)
      div <- 1
      x0 <- 0
      x1 <- -3*tickH
      y0 <- page_height
      y1 <- page_height + 3*tickH
      xsegs <- segmentsGrob(x0 = seq(0, page_width, div), y0 = y0,
                            x1 = seq(0, page_width, div), y1 = y1, default.units = page_units, gp = gpar(col = "black"))
      ysegs <- segmentsGrob(x0 = x0, y0 = seq(0, page_height, div),
                            x1 = x1, y1 = seq(0, page_height, div), default.units = "native", gp = gpar(col = "black"))

      div2 <- 5
      v2 <- segmentsGrob(x0 = seq(0, page_width, div2), y0 = page_height + 3*tickH,
                         x1 = seq(0, page_width, div2), y1 = page_height + 4*tickH, default.units = page_units, gp = gpar(col = "black"))
      h2 <- segmentsGrob(x0 = -3*tickH, y0 = seq(0, page_height, div2),
                         x1 = -4*tickH, y1 = seq(0, page_height, div2), default.units = "native", gp = gpar(col = "black"))
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v2), envir = bbEnv)
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h2), envir = bbEnv)

      div3 <- 10
      v3 <- segmentsGrob(x0 = seq(0, page_width, div3), y0 = page_height + 4*tickH,
                         x1 = seq(0, page_width, div3), y1 = page_height + 5*tickH, default.units = page_units, gp = gpar(col = "black"))
      h3 <- segmentsGrob(x0 = -4*tickH, y0 = seq(0, page_height, div3),
                         x1 = -5*tickH, y1 = seq(0, page_height, div3), default.units = "native", gp = gpar(col = "black"))
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v3), envir = bbEnv)
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h3), envir = bbEnv)

      hLabel <- textGrob(label = seq(0, page_width, div3), x = seq(0, page_width, div3), y = page_height + 4*tickH, vjust = -0.5, default.units = page_units)
      vLabel <- textGrob(label = seq(0, page_height, div3), x = -4*tickH, y = seq(0, page_height, div3), hjust = 1.5, default.units = "native")

    } else {

      tickH <- convertUnit(x = unit(1/32, "inches"), unitTo = page_units, valueOnly = T)
      div <- 1/10
      x0 <- 0
      x1 <- -3*tickH
      y0 <- page_height
      y1 <- page_height + 3*tickH
      xsegs <- segmentsGrob(x0 = seq(0, page_width, div), y0 = y0,
                            x1 = seq(0, page_width, div), y1 = y1, default.units = page_units, gp = gpar(col = "black"))
      ysegs <- segmentsGrob(x0 = x0, y0 = seq(0, page_height, div),
                            x1 = x1, y1 = seq(0, page_height, div), default.units = "native", gp = gpar(col = "black"))


      div2 <- 1/2
      v2 <- segmentsGrob(x0 = seq(0, page_width, div2), y0 = page_height + 3*tickH,
                         x1 = seq(0, page_width, div2), y1 = page_height + 4*tickH, default.units = page_units, gp = gpar(col = "black"))
      h2 <- segmentsGrob(x0 = -3*tickH, y0 = seq(0, page_height, div2),
                         x1 = -4*tickH, y1 = seq(0, page_height, div2), default.units = "native", gp = gpar(col = "black"))
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v2), envir = bbEnv)
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h2), envir = bbEnv)

      div3 <- 1
      v3 <- segmentsGrob(x0 = seq(0, page_width, div3), y0 = page_height + 4*tickH,
                         x1 = seq(0, page_width, div3), y1 = page_height + 5*tickH, default.units = page_units, gp = gpar(col = "black"))
      h3 <- segmentsGrob(x0 = -4*tickH, y0 = seq(0, page_height, div3),
                         x1 = -5*tickH, y1 = seq(0, page_height, div3), default.units = "native", gp = gpar(col = "black"))
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = v3), envir = bbEnv)
      assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = h3), envir = bbEnv)

      hLabel <- textGrob(label = seq(0, page_width, div3), x = seq(0, page_width, div3), y = page_height + 4*tickH, vjust = -0.5, default.units = page_units)
      vLabel <- textGrob(label = seq(0, page_height, div3), x = -4*tickH, y = seq(0, page_height, div3), hjust = 1.5, default.units = "native")

    }

    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = xsegs), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = ysegs), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = hLabel), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = vLabel), envir = bbEnv)


    ## Unit annotation
    unitLabel <- textGrob(label = page_units, x = 0, y = page_height, hjust = 1.75, vjust = -1.5, default.units = page_units, just = c("right", "bottom"))
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = unitLabel), envir = bbEnv)


  }

  # ======================================================================================================================================================================================
  # X AND Y GRIDLINES
  # ======================================================================================================================================================================================

  ## Draw x and y gridlines
  tryCatch(expr = {

    xGrid <- segmentsGrob(x0 = seq(0, page_width, bb_page$xgrid), y0 = 0,
                          x1 = seq(0, page_width, bb_page$xgrid), y1 = bb_page$height,
                          default.units = page_units, gp = gpar(col = "grey50", lty = 2, lwd = 0.5))

    yGrid <- segmentsGrob(x0 = 0, y0 = seq(0, page_height, bb_page$ygrid),
                          x1 = page_width, y1 = seq(0, page_height, bb_page$ygrid),
                          default.units = "native", gp = gpar(col = "grey50", lty = 2, lwd = 0.5))

    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = xGrid), envir = bbEnv)
    assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = yGrid), envir = bbEnv)

  }, error = function(e) return())

  # ======================================================================================================================================================================================
  # DRAW GROBS
  # ======================================================================================================================================================================================

  grid.draw(get("guide_grobs", envir = bbEnv))
  downViewport("bb_page")

}
