#' plots a Manhattan plot
#'
#' @param bedfile bedfile for Manhattan plot
#' @param chrom chromosome of region to be plotted
#' @param chromstart chromstart of region to be plotted
#' @param chromend chromend of region to be plotted
#' @param col single colors, vector of colors, or color palette for coloring points
#' @param ymax fraction of max y value to set as height of plot
#' @param ylim numeric vector of length 2, giving the y coordinates range
#' @param x A unit object specifying x-location
#' @param y A unit object speicifying y-location
#' @param width A unit object specifying width
#' @param height A unit object specifying height
#' @param just justification of plot
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @export
bb_plotManhattan <- function(bedfile, chrom, chromstart, chromend, col, ymax = 1, ylim = NULL, x = NULL, y = NULL, width = NULL, height = NULL, just = "center", draw = TRUE,  ... ){

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

  man_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, x = x, y = y, width = width, height = height, justification = just,
                             colors = col, ylim = ylim, ymax = ymax, grobs = NULL), class = "bb_manhattan")

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

  ## Get viewport name
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_manhattan", length(grep(pattern = "bb_manhattan", x = current_viewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(0.25, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = c(chromstart, chromend), yscale = c(man_plot$ylim[1], man_plot$ylim[2]),
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = man_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(chromstart, chromend), yscale = c(man_plot$ylim[1], man_plot$ylim[2]),
                   just = just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  points <- pointsGrob(x = bedfile[,2], y = -log10(bedfile[,5]), pch = 19, gp = gpar(fill = col), default.units = "native")
  points_grobs <- gTree(vp = vp, children = gList(points))

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(points_grobs)

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  man_plot$grobs <- points_grobs

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(man_plot)

}
