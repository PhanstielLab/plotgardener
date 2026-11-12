#' Export the current plotgardener page without calling dev.off()
#'
#' Grab the current page and write it to a file using a new graphics device.
#' The device is opened and closed internally, so users do not need to call
#' \code{dev.off()}.
#'
#' @param file Output file path. The device is inferred from the extension
#' (`.pdf`, `.png`, `.tiff`, `.jpeg`, or `.svg`).
#' @param width Output width. Defaults to the current page width.
#' @param height Output height. Defaults to the current page height.
#' @param units Units for \code{width} and \code{height}. Defaults to the
#' current page units.
#' @param dpi Resolution for raster devices. Default value is
#' \code{dpi = 300}.
#' @param bg Background color for the output device.
#' Default value is \code{bg = "white"}.
#' @param guides Logical value indicating whether to include page guides in
#' the export. Default value is \code{guides = TRUE}.
#' @param device Optional device name. If \code{NULL}, inferred from
#' \code{file}.
#' @param ... Additional arguments passed to the graphics device.
#'
#' @return Invisibly returns the output file path.
#'
#' @examples
#' \dontrun{
#' pageCreate(width = 3, height = 3, default.units = "inches")
#' plotRect(x = 0.5, y = 0.5, width = 2, height = 2, fill = "green")
#' pageExport("figure.pdf")
#' pageExport("figure.png", dpi = 300)
#' }
#' @export
pageExport <- function(file, width = NULL, height = NULL, units = NULL,
                       dpi = 300, bg = "white", guides = TRUE,
                       device = NULL, ...) {

    if (missing(file) || is.null(file)) {
        stop("argument \"file\" is missing, with no default.", call. = FALSE)
    }
    if (!exists("guide_grobs", envir = pgEnv)) {
        stop("Cannot export without a `plotgardener` page.", call. = FALSE)
    }

    if (is.null(width)) {
        width <- get("page_width", envir = pgEnv)
    }
    if (is.null(height)) {
        height <- get("page_height", envir = pgEnv)
    }
    if (is.null(units)) {
        units <- get("page_units", envir = pgEnv)
    }
    rasterUnits <- if (identical(units, "inches")) "in" else units
    if (!rasterUnits %in% c("in", "cm", "mm")) {
        rasterUnits <- "in"
    }

    ext <- tolower(tools::file_ext(file))
    if (is.null(device)) {
        device <- switch(ext,
            pdf = "pdf",
            png = "png",
            tif = "tiff",
            tiff = "tiff",
            jpg = "jpeg",
            jpeg = "jpeg",
            svg = "svg",
            stop("Unsupported output extension: .", ext, call. = FALSE)
        )
    }

    hideGuides <- !isTRUE(guides) ||
        isTRUE(get("guides_hidden", envir = pgEnv))
    removed <- get("removed_plots", envir = pgEnv)
    needsRemoval <- hideGuides || length(removed) > 0

    ## `wrap.grobs = TRUE` preserves repeated grob names faithfully. When
    ## deferred removals are requested, use the unwrapped tree so individual
    ## guide/plot grobs can be removed by name.
    pageGrob <- grid.grab(wrap.grobs = !needsRemoval)
    if (needsRemoval) {
        if (hideGuides && "guide_grobs" %in% pageGrob$childrenOrder) {
            pageGrob <- removeGrob(pageGrob, "guide_grobs")
        }
        for (grobName in removed) {
            if (grobName %in% pageGrob$childrenOrder) {
                pageGrob <- removeGrob(pageGrob, grobName)
            }
        }
    }

    oldDevice <- grDevices::dev.cur()
    if (device == "pdf") {
        grDevices::pdf(file = file, width = width, height = height, bg = bg, ...)
    } else if (device == "png") {
        grDevices::png(file = file, width = width, height = height,
                       units = rasterUnits, res = dpi, bg = bg, ...)
    } else if (device == "tiff") {
        grDevices::tiff(file = file, width = width, height = height,
                        units = rasterUnits, res = dpi, bg = bg, ...)
    } else if (device == "jpeg") {
        grDevices::jpeg(file = file, width = width, height = height,
                        units = rasterUnits, res = dpi, bg = bg, ...)
    } else if (device == "svg") {
        grDevices::svg(file = file, width = width, height = height, bg = bg, ...)
    } else {
        stop("Unsupported device: ", device, call. = FALSE)
    }
    on.exit({
        if (grDevices::dev.cur() != oldDevice) {
            grDevices::dev.off()
        }
    }, add = TRUE)

    grid.newpage()
    grid.draw(pageGrob)

    invisible(file)
}
