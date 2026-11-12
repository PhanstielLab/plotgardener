#' @keywords internal
.modifyAxisGrob <- function(axisGrob, main = TRUE, tcl = 0.5,
                            label.offset = 0.5,
                            axisLine.lwd = NULL, axisLine.col = NULL,
                            axisLine.lty = NULL,
                            tick.lwd = NULL, tick.col = NULL,
                            tick.lty = NULL,
                            axisLine = TRUE) {
    direction <- if (isTRUE(main)) -1 else 1

    ## Axis line
    if (!axisLine) {
        axisGrob <- removeGrob(axisGrob, "major")
    } else {
        gp <- axisGrob$children$major$gp
        if (!is.null(axisLine.lwd)) gp$lwd <- axisLine.lwd
        if (!is.null(axisLine.col)) gp$col <- axisLine.col
        if (!is.null(axisLine.lty)) gp$lty <- axisLine.lty
        axisGrob$children$major$gp <- gp
    }

    ## Ticks
    if (!is.null(axisGrob$children$ticks)) {
        gp <- axisGrob$children$ticks$gp
        if (!is.null(tick.lwd)) gp$lwd <- tick.lwd
        if (!is.null(tick.col)) gp$col <- tick.col
        if (!is.null(tick.lty)) gp$lty <- tick.lty
        axisGrob$children$ticks$gp <- gp

        if (inherits(axisGrob, "yaxis")) {
            axisGrob$children$ticks$x1 <- unit(direction * tcl, "lines")
        } else {
            axisGrob$children$ticks$y1 <- unit(direction * tcl, "lines")
        }
    }

    ## Labels
    if (!is.null(axisGrob$children$labels)) {
        if (inherits(axisGrob, "yaxis")) {
            axisGrob$children$labels$x <- unit(
                direction * (tcl + label.offset), "lines"
            )
        } else {
            axisGrob$children$labels$y <- unit(
                direction * (tcl + label.offset), "lines"
            )
        }
    }

    axisGrob
}
