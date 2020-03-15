#' draws a horizontal guideline at a specified y-coordinate
#'
#' @param y A numeric or unit object specifying y-coordinate of guide
#' @param col color of guideline
#' @param default.units A string indicating the default units to use if y is only given as numeric vectors

#' @export
bb_hGuide <- function(y, col = "grey55", default.units = "inches"){

  ## Get the names of the current viewports
  # current_viewports <- lapply(current.vpTree()$children, viewport_name)
  #
  # if (!"bb_page" %in% current_viewports){
  #
  #   stop("Must make a BentoBox page with bb_makePage() before adding additional guidelines.")
  #
  # }
  if (class(y) != "unit"){

    if (!is.numeric(y)){

      stop("y-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(default.units)){

      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    y <- unit(y, default.units)
  }

  y <- convertY(y, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)

  guide <- grid.segments(x0 = unit(0, units = "npc"), x1 = unit(1, units = "npc"), y0 = get("page_height", envir = bbEnv) - y, y1 = get("page_height", envir = bbEnv) - y,
                           default.units = get("page_units", envir = bbEnv), gp = gpar(col = col))
  assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = guide), envir = bbEnv)

}
