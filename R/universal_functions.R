## Define a function to turn gTree child into a gPath
convert_gpath <- function(grob){

  ## Get the name of the grob
  name <- grob$name

  ## Turn it into a gPath
  gpath <- gPath(name)

  return(gpath)

}

<<<<<<< HEAD
viewport_name <- function(viewport){
=======
## Define a function that adds grobs to a gTree in bbEnv and draws them
appendGrob <- function(grob, gtree){
  assign(gtree, addGrob(gTree = get(gtree, envir = bbEnv), child = grob), envir = bbEnv)
  grid.draw(grob)
}

>>>>>>> master

  return(viewport$name)

}

## Define a function to convert plot x and y into center of plot based on justification
adjust_coords <- function(plot, page_units, page_height){

  plot_y <- convertY(unit(plot$y, units = plot$units), unitTo = page_units, valueOnly = TRUE)
  plot_y <- page_height - plot_y
  plot_y <- convertY(unit(plot_y, units = page_units), unitTo = plot$units, valueOnly = TRUE)

  if (length(plot$just == 2)){

    if ("left" %in% plot$just & "center" %in% plot$just){
      ## convert the x-coordinate only
      plot_x <- plot$x + (0.5 * plot$width)
    } else if ("right" %in% plot$just & "center" %in% plot$just){
      ## convert the x-coordinate only
      plot_x <- plot$x - (0.5 * plot$width)
    } else if ("center" %in% plot$just & "bottom" %in% plot$just){
      ## convert the y-coordinate only
      plot_x <- plot$x
      plot_y <- plot_y + (0.5 * plot$height)
    } else if ("center" %in% plot$just & "top" %in% plot$just){
      ## convert the y-coordinate only
      plot_x <- plot$x
      plot_y <- plot_y - (0.5 * plot$height)
    } else if ("left" %in% plot$just & "top" %in% plot$just){
      ## convert x-coordinate and y-coordinate
      plot_x <- plot$x + (0.5 * plot$width)
      plot_y <- plot_y - (0.5 * plot$height)
    } else if ("right" %in% plot$just & "top" %in% plot$just){
      ## convert x-coordinate and y-coordinate
      plot_x <- plot$x - (0.5 * plot$width)
      plot_y <- plot_y - (0.5 * plot$height)
    } else if ("left" %in% plot$just & "bottom" %in% plot$just){
      ## convert x-coordinate and y-coordinate
      plot_x <- plot$x + (0.5 * plot$width)
      plot_y <- plot_y + (0.5 * plot$height)
    } else if ("right" %in% plot$just & "bottom" %in% plot$just){
      ## convert x-coordinate and y-coordinate
      plot_x <- plot$x - (0.5 * plot$width)
      plot_y <- plot_y + (0.5 * plot$height)
    } else {
      ## no conversion
      plot_x <- plot$x
    }

  } else if (length(plot$just == 1)){

    if (plot$just == "left"){
      ## convert the x-coordinate only
      plot_x <- plot$x + (0.5 * plot$width)
    } else if (plot$just == "right"){
      ## convert the x-coordinate only
      plot_x <- plot$x - (0.5 * plot$width)
    } else if (plot$just == "bottom"){
      ## convert the y-coordinate only
      plot_x <- plot$x
      plot_y <- plot_y + (0.5 * plot$height)
    } else if (plot$just == "top"){
      ## convert the y-coordinate only
      plot_x <- plot$x
      plot_y <- plot_y - (0.5 * plot$height)
    } else {
      ## no conversion
      plot_x <- plot$x
    }

  }

  return(list(plot_x, plot_y))

}

## Define a function to convert annotation viewport x and y into center of annotation based on justification
adjust_vpCoords <- function(plot, viewport, page_units, page_height){
  #
  # plot_y <- convertY(unit(plot$y, units = plot$units), unitTo = page_units, valueOnly = TRUE)
  # plot_y <- page_height - plot_y
  # plot_y <- convertY(unit(plot_y, units = page_units), unitTo = plot$units, valueOnly = TRUE)

  vp_y <- viewport$y

  if (length(viewport$justification == 2)){

    if ("left" %in% viewport$justification & "center" %in% viewport$justification){

      ## convert the x-coordinate only
      vp_x <- viewport$x + (0.5 * viewport$width)
    } else if ("right" %in% viewport$justification & "center" %in% viewport$justification){
      ## convert the x-coordinate only
      vp_x <- viewport$x - (0.5 * viewport$width)
    } else if ("center" %in% viewport$justification & "bottom" %in% viewport$justification){
      ## convert the y-coordinate only
      vp_x <- viewport$x
      vp_y <- vp_y + (0.5 * viewport$height)
    } else if ("center" %in% viewport$justification & "top" %in% viewport$justification){
      ## convert the y-coordinate only
      vp_x <- viewport$x
      vp_y <- vp_y - (0.5 * viewport$height)
    } else if ("left" %in% viewport$justification & "top" %in% viewport$justification){
      ## convert x-coordinate and y-coordinate
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- vp_y - (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "top" %in% viewport$justification){
      ## convert x-coordinate and y-coordinate
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- vp_y - (0.5 * viewport$height)
    } else if ("left" %in% pviewport$justification & "bottom" %in% viewport$justification){
      ## convert x-coordinate and y-coordinate
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- vp_y + (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "bottom" %in% viewport$justification){
      ## convert x-coordinate and y-coordinate
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- vp_y + (0.5 * viewport$height)
    } else {
      ## no conversion
      vp_x <- viewport$x
    }

  } else if (length(viewport$justification == 1)){

    if (viewport$justification == "left"){
      ## convert the x-coordinate only
      vp_x <- viewport$x + (0.5 * viewport$width)
    } else if (viewport$justification == "right"){
      ## convert the x-coordinate only
      vp_x <- viewport$x - (0.5 * viewport$width)
    } else if (viewport$justification == "bottom"){
      ## convert the y-coordinate only
      vp_x <- viewport$x
      vp_y <- vp_y + (0.5 * viewport$height)
    } else if (viewport$justification == "top"){
      ## convert the y-coordinate only
      vp_x <- viewport$x
      vp_y <- vp_y - (0.5 * viewport$height)
    } else {
      ## no conversion
      vp_x <- viewport$x
    }

  }

  return(list(vp_x, vp_y))

}

# vp_topLeft <- function(viewport, page_units, page_height){
#
#
#
# }



## Define a function to convert to page units
convert_page <- function(object){

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  ## Convert x and y coordinates and height and width to same page_units
  old_x <- unit(object$x, units = object$units)
  old_y <- unit(object$y, units = object$units)
  old_height <- unit(object$height, units = object$units)
  old_width <- unit(object$width, units = object$units)
  new_x <- convertX(old_x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(old_y, unitTo = page_units, valueOnly = TRUE)
  new_height <- convertHeight(old_height, unitTo = page_units, valueOnly = TRUE)
  new_width <- convertWidth(old_width, unitTo = page_units, valueOnly = TRUE)

  object$x <- new_x
  object$y <- new_y
  object$height <- new_height
  object$width <- new_width


  return(list(object, page_height, page_units))

}

## Define a function to make sure a bb_page viewport exists
check_bbpage <- function(){

  ## Get the names of the current viewports
  current_viewports <- lapply(current.vpTree()$children, viewport_name)

  if (!"bb_page" %in% current_viewports){

    stop("Must make a BentoBox page with bb_makePage() before plotting.")

  }

}




