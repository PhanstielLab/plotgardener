#' draws a vertical guideline at a specified x-coordinate
#'
#' @param x A numeric or unit object specifying x-coordinate of guide
#' @param col color of guideline
#' @param default.units A string indicating the default units to use if x is only given as numeric vectors

#' @export
bb_vGuide <- function(x, col = "grey55", default.units = "inches"){

  ## Get the names of the current viewports
  # current_viewports <- lapply(current.vpTree()$children, viewport_name)
  #
  # if (!"bb_page" %in% current_viewports){
  #
  #   stop("Must make a BentoBox page with bb_makePage() before adding additional guidelines.")
  #
  # }
  if (class(x) != "unit"){

    if (!is.numeric(x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    x <- unit(x, default.units)
  }

  guide <- grid.segments(x0 = x, x1 = x, y0 = unit(0, units = "npc"), y1 = unit(1, units = "npc"),
                           gp = gpar(col = col))
  assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = guide), envir = bbEnv)

}
