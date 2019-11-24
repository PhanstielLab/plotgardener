## Define a function to turn gTree child into a gPath
convert_gpath <- function(grob){

  ## Get the name of the grob
  name <- grob$name

  ## Turn it into a gPath
  gpath <- gPath(name)

  return(gpath)

}

## Define a function that adds grobs to a gTree in bbEnv and draws them
appendGrob <- function(grob, gtree){
  assign(gtree, addGrob(gTree = get(gtree, envir = bbEnv), child = grob), envir = bbEnv)
  grid.draw(grob)
}




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






