#' creates a group of grobs
#'
#' @param plots a list of plot objects returned from BentoBox functions
#'
#' @return Function will return a "bb_group" object
#'
#' @export
#'
bb_groupPlots <- function(plots){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function to get the left-most coordinate
  get_Left <- function(plot){

    center_x <- adjust_vpCoords(plot$grobs$vp)[[1]]
    left <- convertX(center_x - (0.5 * plot$grobs$vp$width), unitTo = get("page_units", envir = bbEnv))


    return(left)

  }

  ## Define a function to get the right-most coordinate
  get_Right <- function(plot){

    center_x <- adjust_vpCoords(plot$grobs$vp)[[1]]
    right <- convertX(center_x + (0.5 * plot$grobs$vp$width), unitTo = get("page_units", envir = bbEnv))

    return(right)

  }

  ## Define a function to get the top-most coordinate
  get_Top <- function(plot){

    center_y <- adjust_vpCoords(plot$grobs$vp)[[2]]
    top <- convertY(center_y + (0.5 * plot$grobs$vp$height), unitTo = get("page_units", envir = bbEnv))

    return(top)

  }

  ## Define a function to get the bottom-most coordinate
  get_Bottom <- function(plot){

    center_y <- adjust_vpCoords(plot$grobs$vp)[[2]]
    bottom <- convertY(center_y - (0.5 * plot$grobs$vp$height), unitTo = get("page_units", envir = bbEnv))

    return(bottom)

  }

  ## Define a function to convert plot viewport x and y coordinates to "native" scale
  convert_native <- function(plot){

    plot_vp <- plot$grobs$vp
    plot_vp$x <- convertX(x = plot_vp$x, unitTo = "native")
    plot_vp$y <- convertY(x = plot_vp$y, unitTo = "native")

    plot$grobs$vp <- plot_vp

    return(plot)
  }

  ## Define a function to convert plot viewport widths and heights to "npc" scale
  convert_npc <- function(plot){

    plot_vp <- plot$grobs$vp
    plot_vp$width <- convertWidth(x = plot_vp$width, unitTo = "npc")
    plot_vp$height <- convertHeight(x = plot_vp$height, unitTo = "npc")

    plot$grobs$vp <- plot_vp

    return(plot)
  }

  ## Define a function to add plot gTree to container gTree
  add_gtree <- function(plot){

    gtree <- plot$grobs

    assign("gtree_container", addGrob(gTree = get("gtree_container", envir = bbEnv), child = gtree), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # EXTRACT PLOTS INFO
  # ======================================================================================================================================================================================

  plotLefts <- unlist(lapply(plots, get_Left))
  leftMost <- unit(min(plotLefts), units = get("page_units", envir = bbEnv))
  plotRights <- unlist(lapply(plots, get_Right))
  rightMost <- unit(max(plotRights), units = get("page_units", envir = bbEnv))
  plotTops <- unlist(lapply(plots, get_Top))
  topMost <- unit(max(plotTops), units = get("page_units", envir = bbEnv))
  plotBottoms <- unlist(lapply(plots, get_Bottom))
  bottomMost <- unit(min(plotBottoms), units = get("page_units", envir = bbEnv))

  # ======================================================================================================================================================================================
  # DETERMINE DIMENSIONS AND SCALES OF CONTAINER VIEWPORT
  # ======================================================================================================================================================================================

  vpHeight <- topMost - bottomMost
  vpWidth <- rightMost - leftMost

  xscale1 <- convertX(x = leftMost, unitTo = "native", valueOnly = T)
  xscale2 <- convertX(x = rightMost, unitTo = "native", valueOnly = T)

  yscale1 <- convertY(x = bottomMost, unitTo = "native", valueOnly = T)
  yscale2 <- convertY(x = topMost, unitTo = "native", valueOnly = T)

  # ======================================================================================================================================================================================
  # INTIALIZE OBJECT
  # ======================================================================================================================================================================================

  plot_group <- structure(list(grobs = NULL), class = "bb_group")

  # ======================================================================================================================================================================================
  # MAKE OUTSIDE VIEWPORT TO CONTAIN ALL PLOTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  current_viewports <- current_viewports()
  vp_name <- paste0("bb_group", length(grep(pattern = "bb_group", x = current_viewports)) + 1)

  ## Define viewport
  container_viewport <- viewport(width = vpWidth, height = vpHeight,
                                 x = leftMost, y = topMost, just = c("left", "top"),
                                 xscale = c(xscale1, xscale2), yscale = c(yscale1, yscale2),
                                 name = vp_name)

  # ======================================================================================================================================================================================
  # ADJUST PLOT VIEWPORT X AND Y COORDINATES TO NATIVE
  # ======================================================================================================================================================================================

  plots <- lapply(plots, convert_native)

  # ======================================================================================================================================================================================
  # ADJUST PLOT VIEWPORT WIDTHS AND HEIGHTS TO NPC
  # ======================================================================================================================================================================================

  pushViewport(container_viewport)

  plots <- lapply(plots, convert_npc)

  upViewport()

  # ======================================================================================================================================================================================
  # CREATE GTREE FOR ENTIRE CONTAINER
  # ======================================================================================================================================================================================

  assign("gtree_container", gTree(vp = container_viewport), envir = bbEnv)

  # ======================================================================================================================================================================================
  # ADD PLOT GTREES TO CONTAINER GTREE
  # ======================================================================================================================================================================================

  invisible(lapply(plots, add_gtree))

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  plot_group$grobs <- get("gtree_container", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(plot_group)

}
