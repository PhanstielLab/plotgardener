#' Reshow guides drawn with \code{bb_pageCreate}, \code{bb_pageGuideHorizontal}, and \code{bb_pageGuideVertical}
#'
#'
#' @seealso \link[BentoBox]{bb_pageCreate}, \link[BentoBox]{bb_pageGuideHorizontal}, \link[BentoBox]{bb_pageGuideVertical}
#' @export
bb_pageGuideShow <- function(){

  ## Errors
  ## Get the names of the current viewports
  # current_viewports <- lapply(current.vpTree()$children, viewport_name)
  # if (!"bb_page" %in% current_viewports){
  #
  #   stop("No BentoBox page guides were previously drawn.")
  #
  # }

  if (length(get("guide_grobs", envir = bbEnv)$children) == 0) {

    stop("No BentoBox page guides were previously drawn.")

  }

  grid.draw(get("guide_grobs", envir = bbEnv))

}
