#' zooms and plots region of Hi-C plot
#'
#' @param hic hic plot to zoom inon
#' @param chromstart chromosome start of region to be zoomed
#' @param chromend chromosome end of region to be zoomed
#' @param x x-coordinate of where to place zoom plot based on just
#' @param y y-coordinate of where to place zoom plot based on just
#' @param width width of zoom plot in specified units
#' @param height height of zoom plot in specified units
#' @param units units of height, width, and x and y location of the zoom plot
#' @param just a string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param lty line type of zoom box; options are "solid" (1), "dashed" (2), "dotted" (3), "dotdash" (4), "longdash" (5), or "twodash" (6)
#' @param lwd line width of zoom box
#' @param col line color of zoom box
#' @param lineSide side of zoom box from which to draw lines; options are NULL, "right", "left", "top", or "bottom"
#'


bb_zoomHic <- function(hic, chromstart, chromend, altchromstart = NULL, altchromend = NULL, x, y, width, height, units, just, lty = "solid", lwd = 1, col = "black",
                       lineSide = "right"){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  #errors: hic plot isn't plotted


  get_zoom_grobs <- function(grob, object, hic){

    left_x <- convertUnit(grob$x[1], unitTo = "native", valueOnly = TRUE)
    right_x <- convertUnit(grob$x[4], unitTo = "native", valueOnly = TRUE)

    bottom_y <- convertUnit(grob$y[1], unitTo = "native", valueOnly = TRUE)
    top_y <- convertUnit(grob$y[2], unitTo = "native", valueOnly = TRUE)

    shifted_chromstart <- object$chromstart - hic$additional_parameters$resolution
    shifted_chromend <- object$chromend + hic$additional_parameters$resolution

    if (left_x >= shifted_chromstart & right_x <= shifted_chromend
        & bottom_y >= shifted_chromstart & top_y <= shifted_chromend){

      return(grob)

    }

  }

  copy_grobs <- function(grob, gtree, letter){

    old_name <- grob$name
    new_name <- paste0(old_name, letter)
    grob$name <- new_name

    assign(new_name, grob, envir = bbEnv)

    ## Add to gtree
    assign(gtree, addGrob(gTree = get(gtree, envir = bbEnv), child = grob), envir = bbEnv)

    return(get(new_name, envir = bbEnv))
  }

  plot_zoomLines <- function(lineSide, object){

    ## convert viewport units to center based on justification
    # adjusted_zoomCoords <- adjust_vpCoords(plot = object, viewport = object$viewport[[1]],
    #                                        page_height = get("page_height", envir = bbEnv), page_units = get("page_units", envir = bbEnv))
    #
    # adjusted_zoomCoords <- adjust_vpCoords(plot = object, viewport = object$viewport[[2]],
    #                                        page_height = get("page_height", envir = bbEnv), page_units = get("page_units", envir = bbEnv))

    ## Convert annotation viewport coordinates to top left based on justification (viewport should already be in page_units)
    annot_topLeft <- vp_topLeft(viewport = object$viewport[[2]])

    ## Go into annotation viewport to convert those coordinates to be compatible with the page
    downViewport(name = viewport_name(viewport = object$viewport[[2]]))

    annot_xCenter <- convertX(unit(0.5 * (object$chromend + object$chromstart), unit = "native"), unitTo = get("page_units", envir = bbEnv))
    annot_yCenter <- convertY(unit(0.5 * (object$chromend + object$chromstart), unit = "native"), unitTo = get("page_units", envir = bbEnv))

    upViewport()

    annot_xCenter <- annot_xCenter + annot_topLeft[[1]]
    annot_yCenter <- annot_yCenter
    print(annot_xCenter)
    print(annot_yCenter)


    stop("testing")
    if (lineSide == "right"){

      ## top line goes from top right of vpAnnot to top left of vpZoom
      ## bottom line goes from bottom right of vpAnnot to bottom left of vpZoom

      annot_xRight <- convertUnit(adjusted_annotCoords[[1]] + 0.5 * object$viewport[[2]]$width, unitTo = object$units)
      annot_yTop <- convertUnit(adjusted_annotCoords[[2]] + 0.5 * object$viewport[[2]]$height, unitTo = object$units)
      annot_yBottom <- convertUnit(adjusted_annotCoords[[2]] - 0.5 * object$viewport[[2]]$height, unitTo = object$units)


      print(annot_xRight)
      print(annot_yTop)
      print(annot_yBottom)

      zoom_xLeft <- convertUnit(adjusted_zoomCoords[[1]] - 0.5 * object$viewport[[2]]$width, unit = object$units)
      zoom_yTop <- convertUnit(adjusted_zoomCoords[[2]] + 0.5 * object$viewport[[2]]$height, unit = object$units)
      zoom_yBottom <- convertUnit(adjusted_zoomCoords[[2]] - 0.5 * object$viewport[[2]]$height, unit = object$units)

      print(zoom_xLeft)
      print(zoom_yTop)
      print(zoom_yBottom)

      topLine <- segmentsGrob(x0 = annot_xRight, y0 = annot_yTop, x1 = zoom_xLeft, y1 = zoom_yTop,
                              gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col))
      bottomLine <- segmentsGrob(x0 = annot_xRight, y0 = annot_yBottom, x1 = zoom_xLeft, y1 = zoom_yBottom,
                                 gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col))

    } else if (lineSide == "left"){

      ## top line goes from top left of vpAnnot to top right of vpZoom
      ## bottom line goes from bottom left of vpAnnot to bottom right of vpZoom

      annot_xLeft <- adjusted_annotCoords[[1]] - 0.5 * object$viewport[[2]]$width
      annot_yTop <- adjusted_annotCoords[[2]] + 0.5 * object$viewport[[2]]$height
      annot_yBottom <- adjusted_annotCoords[[2]] - 0.5 * object$viewport[[2]]$height

      zoom_xRight <- adjusted_zoomCoords[[1]] + 0.5 * object$viewport[[2]]$width
      zoom_yTop <- adjusted_zoomCoords[[2]] + 0.5 * object$viewport[[2]]$height
      zoom_yBottom <- adjusted_zoomCoords[[2]] - 0.5 * object$viewport[[2]]$height

      topLine <- segmentsGrob(x0 = annot_xLeft, y0 = annot_yTop, x1 = zoom_xRight, y1 = zoom_yTop,
                              gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col))
      bottomLine <- segmentsGrob(x0 = annot_xLeft, y0 = annot_yBottom, x1 = zoom_xRight, y1 = zoom_yBottom,
                                 gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col))

    } else if (lineSide == "top"){

      ## top line goes from top left of vpAnnot to bottom left of vpZoom
      ## bottom line goes from top right of vpAnnot to bottom right of vpZoom

      annot_xRight <- adjusted_annotCoords[[1]] + 0.5 * object$viewport[[2]]$width
      annot_xLeft <- adjusted_annotCoords[[1]] - 0.5 * object$viewport[[2]]$width
      annot_yTop <- adjusted_annotCoords[[2]] + 0.5 * object$viewport[[2]]$height

      zoom_xRight <- adjusted_zoomCoords[[1]] + 0.5 * object$viewport[[2]]$width
      zoom_xLeft <- adjusted_zoomCoords[[1]] - 0.5 * object$viewport[[2]]$width
      zoom_yBottom <- adjusted_zoomCoords[[2]] - 0.5 * object$viewport[[2]]$height

      topLine <- segmentsGrob(x0 = annot_xLeft, y0 = annot_yTop, x1 = zoom_xLeft, y1 = zoom_yBottom,
                              gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col))
      bottomLine <- segmentsGrob(x0 = annot_xRight, y0 = annot_yTop, x1 = zoom_xRight, y1 = zoom_yBottom,
                                 gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col))

    } else if (lineSide == "bottom"){

      ## top line goes from bottom left of vpAnnot to top left of vpZoom
      ## bottom line goes from bottom right of vpAnnot to top right of vpZoom

      annot_xRight <- adjusted_annotCoords[[1]] + 0.5 * object$viewport[[2]]$width
      annot_xLeft <- adjusted_annotCoords[[1]] - 0.5 * object$viewport[[2]]$width
      annot_yBottom <- adjusted_annotCoords[[2]] - 0.5 * object$viewport[[2]]$height

      zoom_xRight <- adjusted_zoomCoords[[1]] + 0.5 * object$viewport[[2]]$width
      zoom_xLeft <- adjusted_zoomCoords[[1]] - 0.5 * object$viewport[[2]]$width
      zoom_yTop <- adjusted_zoomCoords[[2]] + 0.5 * object$viewport[[2]]$height

      topLine <- segmentsGrob(x0 = annot_xLeft, y0 = annot_yBottom, x1 = zoom_xLeft, y1 = zoom_yTop,
                              gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col))
      bottomLine <- segmentsGrob(x0 = annot_xRight, y0 = annot_yBottom, x1 = zoom_xRight, y1 = zoom_yTop,
                                 gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col))

    }

    grid.draw(topLine)
    grid.draw(bottomLine)
    assign("zoom_grobs", addGrob(gTree = get("zoom_grobs", envir = bbEnv), child = topLine), envir = bbEnv)
    assign("zoom_grobs", addGrob(gTree = get("zoom_grobs", envir = bbEnv), child = bottomLine), envir = bbEnv)


  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  hic_zoom <- structure(list(chrom = hic$chrom, chromstart = chromstart, chromend = chromend, altchrom = NULL, altchromstart = NULL, altchromend = NULL, x = x, y = y, width = width,
                               height = height, units = units, justification = just, grobs = NULL, viewport = NULL,
                               gpar = list(lty = lty, lwd = lwd, col = col)), class = "bb_hicZoom")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage()

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE OF GROBS
  # ======================================================================================================================================================================================

  assign("zoom_grobs", gTree(name = "zoom_grobs"), envir = bbEnv)

  # ======================================================================================================================================================================================
  # GET GROBS IN SPECIFIED REGION
  # ======================================================================================================================================================================================

  ## Get whether grobs are within region
  zoom_grobs <- lapply(hic$grobs, get_zoom_grobs, object = hic_zoom, hic = hic)

  ## Remove null values
  zoom_grobs <- zoom_grobs[!unlist(lapply(zoom_grobs, is.null))]

  ## Change grob names/assignments so they aren't the same as the hic plot itself/other zoom plots and add to gtree
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  zoom_grobs <- lapply(zoom_grobs, copy_grobs, gtree = "zoom_grobs", letter = letters[length(grep(pattern = "bb_hicZoom", x = current_viewports)) + 1])
  assign("ZOOM_GROBS", zoom_grobs, envir = globalenv())

  # ======================================================================================================================================================================================
  # ZOOM PLOT VIEWPORT
  # ======================================================================================================================================================================================

  ## Make name of viewport
  vp_name1 <- paste0("bb_hicZoom", length(grep(pattern = "bb_hicZoom", x = current_viewports)) + 1)

  ## Get viewport xscale and yscale

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = hic_zoom)

  ## Make viewport based on hic input viewport
  vp1 <- viewport(height = unit(page_coords[[1]]$height, page_coords[[3]]), width = unit(page_coords[[1]]$width, page_coords[[3]]),
                 x = unit(page_coords[[1]]$x, page_coords[[3]]), y = unit((page_coords[[2]]-page_coords[[1]]$y), page_coords[[3]]),
                 clip = "on", xscale = c(chromstart, chromend), yscale = c(chromstart, chromend), just = just, name = vp_name1)

  pushViewport(vp1)


  # ======================================================================================================================================================================================
  # PLOT ZOOM
  # ======================================================================================================================================================================================

  invisible(lapply(zoom_grobs, grid.draw))

  ## Go back to root viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ANNOTATION VIEWPORT
  # ======================================================================================================================================================================================

  vp_name2 <- paste0("bb_zoomAnnotation", length(grep(pattern = "bb_zoomAnnotation", x = current_viewports)) + 1)

  ## Make viewport based on hic input viewport
  vp2 <- viewport(height = hic$viewport$height, width = hic$viewport$width,
                 x = hic$viewport$x, y = hic$viewport$y, clip = "on", xscale = hic$viewport$xscale, yscale = hic$viewport$yscale, just = hic$viewport$justification,
                 name = vp_name2)
  pushViewport(vp2)

  hic_zoom$viewport <- vpList(vp1, vp2)

  # ======================================================================================================================================================================================
  # PLOT ANNOTATION
  # ======================================================================================================================================================================================

  ## Draw box around region
  annot <- rectGrob(x = chromstart, y = chromstart, just = c("left", "bottom"), width = (chromend - chromstart), height = (chromend - chromstart),
           default.units = "native", gp = gpar(fill = NA))
  grid.draw(annot)

  ## Add to gtree
  assign("zoom_grobs", addGrob(gTree = get("zoom_grobs", envir = bbEnv), child = annot), envir = bbEnv)

  ## Go back to root viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ZOOM LINES
  # ======================================================================================================================================================================================
  if (!is.null(lineSide)){

    plot_zoomLines(lineSide = lineSide, object = hic_zoom)

  }


  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  hic_zoom$grobs <- get("zoom_grobs", envir = bbEnv)$children

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(hic_zoom)

}
