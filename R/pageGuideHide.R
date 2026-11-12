#' Remove guides from a plotgardener page
#'
#' Hides page guides without calling `grid.remove()` on the active device.
#' The hidden state is applied by \code{pageExport()} and \code{pageRedraw()},
#' so exporting a PDF or PNG does not create a second page.
#'
#' @param redraw Logical value indicating whether to redraw immediately on a
#' screen device. File devices defer the change to \code{pageExport()}.
#' Default value is \code{redraw = TRUE}.
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
#' ## Export with guides omitted
#' \dontrun{
#' pageExport("figure.pdf")
#' }
#' @export
pageGuideHide <- function(redraw = TRUE) {
    if (length(get("guide_grobs", envir = pgEnv)$children) == 0) {
        stop("No `plotgardener` page guides exist.")
    }

    assign("guides_hidden", TRUE, envir = pgEnv)
    if (isTRUE(redraw) && !.isFileDevice()) {
        pageRedraw()
    }
    invisible(NULL)
}
