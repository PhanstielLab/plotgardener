#' @keywords internal
.isFileDevice <- function() {
    device <- names(grDevices::dev.cur())
    any(grepl("pdf|png|jpeg|jpg|tiff|svg|postscript",
              device, ignore.case = TRUE))
}

#' Redraw the current plotgardener page
#'
#' Redraw the current page after deferred guide hiding or plot removal. This
#' is primarily useful on screen devices; for file devices, use
#' \code{pageExport()}, which applies the deferred removals to a single page.
#'
#' @return Invisibly returns the grabbed page grob.
#' @export
pageRedraw <- function() {
    if (!exists("guide_grobs", envir = pgEnv)) {
        stop("Cannot redraw without a `plotgardener` page.", call. = FALSE)
    }

    pageGrob <- grid.grab(wrap.grobs = FALSE)
    if (isTRUE(get("guides_hidden", envir = pgEnv)) &&
        "guide_grobs" %in% pageGrob$childrenOrder) {
        pageGrob <- removeGrob(pageGrob, "guide_grobs")
    }
    for (grobName in get("removed_plots", envir = pgEnv)) {
        if (grobName %in% pageGrob$childrenOrder) {
            pageGrob <- removeGrob(pageGrob, grobName)
        }
    }

    grid.newpage()
    grid.draw(pageGrob)
    invisible(pageGrob)
}
