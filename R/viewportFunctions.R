## Define a function to grab the name of a viewport
viewport_name <- function(viewport){

  return(viewport$name)

}

## Define a function to add a plot viewport to internal list of valid "last plotted" viewports
add_bbViewport <- function(vpName){

  bb_vpTree <- get("bb_vpTree", envir = bbEnv)
  bb_vpTree[length(bb_vpTree) + 1] <- vpName
  assign("bb_vpTree", bb_vpTree, envir = bbEnv)

}

## Define a function to get a list of current viewports
current_viewports <- function(){

  if (!"bb_page" %in% names(lapply(current.vpTree()$children, viewport_name))){

    current <- as.list(names(lapply(current.vpTree()$children, viewport_name)))

  } else {

    page_children <- names(lapply(current.vpTree()$children$bb_page$children, viewport_name))
    current <- as.list(page_children)

  }

  return(current)
}

## Define a function to convert viewport x and y into center based on justification
adjust_vpCoords <- function(viewport){
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
    } else if ("left" %in% viewport$justification & "bottom" %in% viewport$justification){
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

## Define a function to change viewport x and y-coordinates to top left based on justification
vp_topLeft <- function(viewport){

  if (length(viewport$justification == 2)){

    if ("left" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y + (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y + (0.5 * viewport$height)
    } else if ("center" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y + (viewport$height)
    } else if ("center" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if ("left" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y
    } else if ("right" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y
    } else if ("left" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y + (viewport$height)
    } else if ("right" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y + (viewport$height)
    } else {
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y + (0.5 * viewport$height)
    }

  } else if (length(viewport$justification == 1)){

    if (viewport$justification == "left"){
      vp_x <- viewport$x
      vp_y <- viewport$y + (0.5 * viewport$height)
    } else if (viewport$justification == "right"){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y + (0.5 * viewport$height)
    } else if (viewport$justification == "bottom"){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y + (viewport$height)
    } else if (viewport$justification == "top"){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y
    } else {
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y + (0.5 * viewport$height)
    }

  }

  return(list(vp_x, vp_y))

}

## Define a function to change viewport x and y-coordinates to bottom left based on justification
vp_bottomLeft <- function(viewport){

  if (length(viewport$justification == 2)){

    if ("left" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if ("center" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if ("center" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else if ("left" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y - (viewport$height)
    } else if ("right" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else if ("left" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y
    } else if ("right" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y
    } else {
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    }

  } else if (length(viewport$justification == 1)){

    if (viewport$justification == "left"){
      vp_x <- viewport$x
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (viewport$justification == "right"){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (viewport$justification == "bottom"){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if (viewport$justification == "top"){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else {
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    }

  }

  return(list(vp_x, vp_y))

}

## Define a function to change viewport x and y-coordinates to bottom right based on justification
vp_bottomRight <- function(viewport){

  if (length(viewport$justification == 2)){

    if ("left" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x + viewport$width
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if ("center" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if ("center" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else if ("left" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x + viewport$width
      vp_y <- viewport$y - (viewport$height)
    } else if ("right" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y - (viewport$height)
    } else if ("left" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x + viewport$width
      vp_y <- viewport$y
    } else if ("right" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y
    } else {
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    }

  } else if (length(viewport$justification == 1)){

    if (viewport$justification == "left"){
      vp_x <- viewport$x + viewport$width
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (viewport$justification == "right"){
      vp_x <- viewport$x
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (viewport$justification == "bottom"){
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if (viewport$justification == "top"){
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else {
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    }

  }

  return(list(vp_x, vp_y))

}

## Define a function to convert to page units
convert_page <- function(object){

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  ## Convert x and y coordinates and height and width to same page_units
  old_x <- object$x
  old_y <- object$y
  old_height <- object$height
  old_width <- object$width
  new_x <- convertX(old_x, unitTo = page_units)
  new_y <- convertY(unit(page_height, units = page_units) - convertY(old_y, unitTo = page_units), unitTo = page_units)
  new_height <- convertHeight(old_height, unitTo = page_units)
  new_width <- convertWidth(old_width, unitTo = page_units)

  object$x <- new_x
  object$y <- new_y
  object$height <- new_height
  object$width <- new_width


  return(object)

}

## Parse y-coordinates for plots relative to bottom of last plot
plot_belowY <- function(y_coord){

  prevPlots <- unlist(get("bb_vpTree", envir = bbEnv))

  if (length(prevPlots) == 0){
    stop("No previous plot detected. Cannot define a \'below\' y-coordinate.", call. = FALSE)
  }

  lastPlot <- prevPlots[length(prevPlots)]
  ## Go into viewport to get coordinates/dimensions
  seekViewport(name = lastPlot)
  x <- current.viewport()$x
  y <- current.viewport()$y
  width <- current.viewport()$width
  height <- current.viewport()$height
  just <- current.viewport()$just

  ## Exit viewport
  upViewport()

  ## Convert to bottom
  bottomCoords <- vp_bottomLeft(viewport(x = x, y = y, width = width, height = height, just = just))

  ## Parse number and convert to page units
  space <- as.numeric(gsub("b", "", y_coord))
  space <- unit(space, get("page_units", envir = bbEnv))

  new_y <- (unit(get("page_height", envir = bbEnv), get("page_units", envir = bbEnv)) - bottomCoords[[2]]) + space
  new_y <- convertY(new_y, unitTo = get("page_units", envir = bbEnv))

  return(new_y)
}
