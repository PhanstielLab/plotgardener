#' places a plot that has been previously created but not drawn
#'
#' @param plot plot to be placed
#' @param x A unit object specifying x-location.
#' @param y A unit object specifying y-location.
#' @param width A unit object specifying width.
#' @param height A unit object specifying height.
#' @param just A string or numeric vector specifying the justification of the plot relative to its (x, y) location
#'
#' @return Function will plot a HiC interaction matrix and return a bb_hicPlot object
#'
#'
#' @export
bb_placePlot <- function(plot, x, y, width, height, just = c("left", "top")){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================
  ## Errors: a bb_page doesn't exist


  ## Define a function that renames grobs and copies to a new gtree
  copy_grobs <- function(grob){

    new_name <- grobName(grob)

    grob$name <- new_name
    assign("new_gtree", addGrob(gTree = get("new_gtree", envir = bbEnv), child = grob), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE PLOT OBJECT
  # ======================================================================================================================================================================================

  object <- plot

  # ======================================================================================================================================================================================
  # UPDATE DIMENSIONS AND COORDINATES OF PLOT OBJECT
  # ======================================================================================================================================================================================

  object$x <- x
  object$y <- y
  object$width <- width
  object$height <- height
  object$justification <- just

  # ======================================================================================================================================================================================
  # DEFINE A NEW VIEWPORT
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = object)

  ## Get viewport name
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0(gsub(pattern = "[0-9]", replacement = "", x = object$grobs$vp$name), length(grep(pattern = gsub(pattern = "[0-9]", replacement = "", x = object$grobs$vp$name), x = current_viewports)) + 1)

  ## Make viewport
  new_vp <- viewport(height = page_coords$height, width = page_coords$width,
                 x = page_coords$x, y = page_coords$y,
                 clip = "on",
                 xscale = object$grobs$vp$xscale, yscale = object$grobs$vp$yscale,
                 just = just,
                 name = vp_name)

  # ======================================================================================================================================================================================
  # RENAME GROBS
  # ======================================================================================================================================================================================

  assign("new_gtree", gTree(vp = new_vp), envir = bbEnv)
  invisible(lapply(object$grobs$children, copy_grobs))

  # ======================================================================================================================================================================================
  # ASSIGN NEW GROBS TO OBJECT
  # ======================================================================================================================================================================================

  object$grobs <- get("new_gtree", envir = bbEnv)

  # ======================================================================================================================================================================================
  # DRAW GROBS
  # ======================================================================================================================================================================================

  grid.draw(object$grobs)

  # ======================================================================================================================================================================================
  # RETURN UPDATED OBJECT
  # ======================================================================================================================================================================================

  return(object)


}
