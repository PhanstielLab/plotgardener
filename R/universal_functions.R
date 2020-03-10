## Define a function to turn gTree child into a gPath
convert_gpath <- function(grob){

  ## Get the name of the grob
  name <- grob$name

  ## Turn it into a gPath
  gpath <- gPath(name)

  return(gpath)

}

## Define a function to grab the name of a viewport
viewport_name <- function(viewport){

  return(viewport$name)

}

## Define a function that gets the children of a group viewport
vp_children <- function(group_name){

  children <- unlist(current.vpTree()$children$bb_page$children[group_name], recursive = F)[[2]]
  return(children)

}

## Define a function to get a list of current viewports
current_viewports <- function(){

  if (!"bb_page" %in% names(lapply(current.vpTree()$children, viewport_name))){

    current <- as.list(names(lapply(current.vpTree()$children, viewport_name)))

  } else {

    ## Check for groups
    page_children <- names(lapply(current.vpTree()$children$bb_page$children, viewport_name))

    if (length(grep(pattern = "bb_group", x = page_children)) > 0){

      group_vps <- as.list(page_children[grep(pattern = "bb_group", x = page_children)])

      group_children <- unlist(lapply(group_vps, vp_children), recursive = F)

      children_vps <- lapply(group_children, viewport_name)

      current <- c(page_children, children_vps)

    } else {

      current <- as.list(page_children)

    }

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
    } else if ("left" %in% pviewport$justification & "bottom" %in% viewport$justification){
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
    } else if ("left" %in% pviewport$justification & "bottom" %in% viewport$justification){
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
  new_y <- unit(page_height, units = page_units) - convertY(old_y, unitTo = page_units)
  new_height <- convertHeight(old_height, unitTo = page_units)
  new_width <- convertWidth(old_width, unitTo = page_units)

  object$x <- new_x
  object$y <- new_y
  object$height <- new_height
  object$width <- new_width


  return(object)

}

## Define a function to make sure a bb_page viewport exists
check_bbpage <- function(error){

  if (!"bb_page" %in% current.vpPath()){

    stop(error, call. = FALSE)

  }

}

## Define a function to check dimensions/placing coordinates
check_placement <- function(object){

  if (attributes(object)$plotted == T){

    ## If giving placement coordinates
    if (!is.null(object$x) | !is.null(object$y)){

      ## 1. Need both an x and y coordinate
      if (!is.null(object$x) & is.null(object$y)){

        stop("Placement detected with y value missing.", call. = FALSE)

      }

      if (!is.null(object$y) & is.null(object$x)){

        stop("Placement detected with x value missing.", call. = FALSE)

      }

      ## 2. Need plot dimensions
      if (is.null(object$width)){

        stop("Placement detected with plot width missing.", call. = FALSE)

      }

      if (is.null(object$height)){

        stop("Placement detected with plot height missing.", call. = FALSE)

      }

      ## 3. Need a bb_page
      check_bbpage(error = "Must make a BentoBox page with bb_makePage() before placing a plot.")

    }

  }

}


check_group_placement <- function(object){

  if (attributes(object)$plotted == T){

    ## If plotting a group, need to provide placement coordinates
    if(is.null(object$x) | is.null(object$y)){

      stop("If placing a group, need to specify x and y coordinates.")

    }

    ## Need a bb_page
    check_bbpage(error = "Must make a BentoBox page with bb_makePage() before placing a plot or group.")

  }

}






