## Define a function to convert based on top of page
convert_coordinates <- function(height, width, x, y, pageheight){

  #given x and y coordinates in the top left
  #y coordinate based as if going from  of page

  ybottom = pageheight - y

  x1 = x + (0.5*width)
  y1 = ybottom - (0.5*height)

  return(list(x1, y1))

}


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

## Define a function that adjusts coordinates back to center based on justification
adjust_just <- function(x, y, height, width, just){

  ## options for just are left, right, centre, center, bottom, top; also in various combinations
  if (just == "left" | just == c("left", "center") | just == c("left", "centre")){
    adjX <- x + (0.5 * width)
    adjY <- y
  }

  if (just == "right" | just == c("right", "center") | just == c("right", "centre")){
    adjX <- x - (0.5 * width)
    adjY <- y
  }

  if (just == "bottom" | just == c("center", "bottom") | just == c("centre", "bottom")){
    adjX <- x
    adjY <- y + (0.5 * height)
  }

  if (just == "top" | just == c("center", "top") | just == c("centre", "top")){
    adjX <- x
    adjY <- y - (0.5 * height)
  }

  if (just == c("left", "top")){
    adjX <- x + (0.5 * width)
    adjY <- y - (0.5 * height)
  }

  if (just == c("right", "top")){
    adjX <- x - (0.5 * width)
    adjY <- y - (0.5 * height)
  }

  if (just == c("left", "bottom")){
    adjX <- x + (0.5 * width)
    adjY <- y + (0.5 * height)
  }

  if (just == c("right", "bottom")){
    adjX <- x - (0.5 * width)
    adjY <- y + (0.5 * height)
  }

  if (just == "center" | just == "centre"| just == c("center", "center") | just == c("centre", "centre") |
      just == c("center", "centre") | just == c("centre", "center")){
    adjX <- x
    adjY <- y
  }

  return(list(adjX, adjY))

}




