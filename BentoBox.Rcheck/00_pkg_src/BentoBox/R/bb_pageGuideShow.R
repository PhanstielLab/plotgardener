#' Reshow guides drawn with \code{bb_pageCreate},
#' \code{bb_pageGuideHorizontal}, and \code{bb_pageGuideVertical}
#'
#' @seealso \link[BentoBox]{bb_pageCreate},
#' \link[BentoBox]{bb_pageGuideHorizontal},
#' \link[BentoBox]{bb_pageGuideVertical}
#'
#' @return None.
#' @export
bb_pageGuideShow <- function() {
    if (length(get("guide_grobs", envir = bbEnv)$children) == 0) {
        stop("No BentoBox page guides were previously drawn.")
    }

    grid.draw(get("guide_grobs", envir = bbEnv))
}
