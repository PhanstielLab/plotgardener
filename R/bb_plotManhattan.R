#' plots a Manhattan plot
#'
#' @param bedfile bedfile for Manhattan plot
#' @param chrom chromosome of region to be plotted
#' @param chromstart chromstart of region to be plotted
#' @param chromend chromend of region to be plotted
#' @param col single colors, vector of colors, or color palette for coloring points
#' @param ymax fraction of max y value to set as height of plot
#' @param ylim numeric vector of length 2, giving the y coordinates range
#' @param x x-coordinate
#' @param y y-coordinate
#' @param width width of plot
#' @param height height of plot
#' @param just justification of plot
#' @param units units of width, height, x, and y-coordinates
#'
#' @export
bb_plotManhattan <- function(bedfile, chrom, chromstart, chromend, col, ymax = 1, ylim = NULL, x, y, width, height, just, units,  ... ){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================
  ## Define a function that calls errors for bb_plotManhattan

  ## Define a function that adjusts the range
  manhattan_range <- function(bedfile, object){

    if (is.null(object$ylim)){

      object$ylim <- c(0, object$ymax * max(-log10(bedfile[,5])))

    }

    return(object)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  man_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, x = x, y = y, width = width, height = height, units = units, justification = just,
                             colors = col, ylim = ylim, ymax = ymax, grobs = NULL, viewport = NULL), class = "bb_manhattan")

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  bedfile <- bedfile[which( bedfile[,1] == chrom),]

  # ======================================================================================================================================================================================
  # Y-LIMITS
  # ======================================================================================================================================================================================

  man_plot <- manhattan_range(bedfile = bedfile, object = man_plot)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = man_plot)

  ## Make viewport name
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_manhattan", length(grep(pattern = "bb_manhattan", x = current_viewports)) + 1)

  ## Define viewport
  vp <- viewport(width = unit(page_coords[[1]]$width, units = page_coords[[3]]), height = unit(page_coords[[1]]$height, units = page_coords[[3]]),
                 x = unit(page_coords[[1]]$x, units = page_coords[[3]]), y = unit((page_coords[[2]]-page_coords[[1]]$y), page_coords[[3]]), xscale = c(chromstart, chromend),
                 yscale = c(man_plot$ylim[1], man_plot$ylim[2]), just = just, name = vp_name, clip = "on")

  man_plot$viewport <- vp
  pushViewport(vp)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  points <- grid.points(x = bedfile[,2], y = -log10(bedfile[,5]), pch = 19, gp = gpar(fill = col), default.units = "native")
  points_grobs <- gTree(name = "points_grobs", children = gList(points))

  ## Go back to root viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  man_plot$grobs <- points_grobs$children

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(man_plot)

}
