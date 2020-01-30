#' adds a plot created with ggplot2 to a BentoBox layout
#'
#' @param plot ggplot object
#' @param x A unit object specifying x-location.
#' @param y A unit object specifying y-location.
#' @param width A unit object specifying width.
#' @param height A unit object specifying height.
#' @param just A string or numeric vector specifying the justification of the plot relative to its (x, y) location
#'
#' @export
bb_plotGG <- function(plot, x, y, width, height, just = c("left", "top")){

  # ======================================================================================================================================================================================
  # INITIALIZE PLOT OBJECT
  # ======================================================================================================================================================================================

  gg_plot <- structure(list(width = width, height = height, x = x, y = y, justification = just), class = "bb_gg")

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = gg_plot)

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_gg", length(grep(pattern = "bb_gg", x = currentViewports)) + 1)

  ## Make viewport for gene track
  vp <- viewport(height = page_coords$height, width = page_coords$width,
                 x = page_coords$x, y = page_coords$y,
                 just = just, name = vp_name)

  # ======================================================================================================================================================================================
  # PRINT GGPLOT
  # ======================================================================================================================================================================================

  print(plot, vp = vp)


  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(gg_plot)

}
