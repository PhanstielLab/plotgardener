#' generates new page for plotting
#'
#'
#' @param height height of page
#' @param width width of page
#' @param units units for height and width; options are same as grid options: "inches", "cm", "mm", and "npc"
#' @param showOutline if TRUE, will show limits of defined page
#'
#' @export
#'
#'
#'
bb_makepage <- function(width = 8.5, height = 11, units = "inches", showOutline = FALSE){

  grid.newpage()

  page_vp <- viewport(x = unit(0.5, "npc"), y = unit(0.5, "npc"), width = unit(width, units = units), height = unit(height, units = units), name = "bb_page")

  pushViewport(page_vp)

  if (showOutline == TRUE){
    grid.rect()
  }
  assign("page_height", height, envir = bbEnv)
  assign("page_units", units, envir = bbEnv)

}
