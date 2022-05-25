#' Remove guides from a plotgardener page
#'
#' @return None.
#'
#' @examples
#' ## Make a page
#' pageCreate(width = 7, height = 4, default.units = "inches")
#'
#' ## Hide page guides
#' pageGuideHide()
#' 
#' #' @note 
#' Please note that due to the implementation of `grid` removal functions,
#' using `pageGuideHide` within a `pdf` call will result in the rendering of a 
#' separate, new page with the plot guides removed. To avoid this artifact,
#' hide guides in the `pageCreate` function call with `showGuides = FALSE`.
#' 
#' @export
pageGuideHide <- function() {
    if (length(get("guide_grobs", envir = pgEnv)$children) == 0) {
        stop("No `plotgardener` page guides exist.")
    }

    ## Get the list of grobs from the guide grobs gTree and convert to gPaths
    grob_list <- lapply(
        get("guide_grobs", envir = pgEnv)$children,
        convert_gpath
    )
    ## Get the last grob
    last_grob <- grob_list[length(grob_list)]

    ## Remove the last grob from the list of grobs
    grob_list <- grob_list[-length(grob_list)]

    ## Remove all grobs except last one, not redrawing each time (to save time)
    invisible(lapply(grob_list, grid.remove, redraw = FALSE))

    ## Remove last grob with redrawing, now removing all the grobs
    invisible(lapply(last_grob, grid.remove))
}
