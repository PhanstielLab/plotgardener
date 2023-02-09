#' Create a page for a plotgardener layout
#' 
#' @usage pageCreate(
#'     width = 8.5,
#'     height = 11,
#'     default.units = "inches",
#'     bg = NA,
#'     xgrid = 0.5,
#'     ygrid = 0.5,
#'     showGuides = TRUE,
#'     params = NULL
#' )
#'
#' @param width A numeric or unit object specifying page width.
#' Default value is \code{width = 8}.
#' @param height A numeric or unit object specifying page height.
#' Default value is \code{height = 11}.
#' @param default.units A string indicating the default units to use if
#' \code{width} or \code{height} are only given as numerics.
#' Default value is \code{default.units = "inches"}.
#' @param bg Character value indicating page background color.
#' Default value is \code{bg = NA}.
#' @param xgrid A numeric indicating the increment by which to place
#' vertical gridlines. Default value is \code{xgrid = 0.5}.
#' @param ygrid A numeric indicating the increment by which to place
#' horizontal gridlines. Default value is \code{ygrid = 0.5}.
#' @param showGuides A logical value indicating whether to draw a
#' black border around the entire page and guiding rulers along the top and
#' left side of the page. Default value is \code{showOutline = TRUE}.
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#'
#' @examples
#' ## Create a 6-inch wide, 4.5-inch high page
#' pageCreate(width = 6, height = 4.5, default.units = "inches")
#'
#' ## Create a 14-cm wide, 10-cm high page
#' pageCreate(width = 14, height = 10, default.units = "cm")
#' @details \code{width} and \code{height} must be specified in the same units.
#'
#' @return None.
#'
#' @export
pageCreate <- function(width = 8.5, height = 11, default.units = "inches",
                        bg = NA, xgrid = 0.5, ygrid = 0.5, showGuides = TRUE,
                        params = NULL) {

    # =========================================================================
    # FUNCTION
    # =========================================================================

    errorcheck_pageCreate <- function(object) {

        ## Width and height can't be negative
        if (as.numeric(object$width) <= 0) {
            stop("`width` cannot be 0 or negative.", call. = FALSE)
        }

        if (as.numeric(object$height) <= 0) {
            stop("`height` cannot be 0 or negative.", call. = FALSE)
        }


        ## Width and height need to be of the same units
        widthUnit <- gsub("[0-9]|[.]", "", object$width)
        heightUnit <- gsub("[0-9]|[.]", "", object$height)
        if (widthUnit != heightUnit) {
            stop("`width` and `height` must be in the same units. ",
                "`width` detected as", "`", widthUnit, "`",
                "and `height` detected as", "`", heightUnit, "`", ".",
                call. = FALSE
            )
        }

        ## Page units must be a valid option
        page_units <- gsub("[0-9]|[.]", "", page$width)
        if (!page_units %in% validUnits) {
            message <- paste(c(
                "Invalid page units. Options for page units are:",
                validUnits
            ),
            collapse = " "
            )
            stop(message, call. = FALSE)
        }
    }

    # =========================================================================
    # MAKE NEW PAGE
    # =========================================================================

    grid.newpage()
    assign("pg_vpTree", list(), envir = pgEnv)

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    page <- parseParams(params = params, 
                        defaultArgs = formals(eval(match.call()[[1]])),
                        declaredArgs = lapply(match.call()[-1], eval.parent, 
                                            n = 2),
                        class = "page")
    
    # =========================================================================
    # DEFAULT UNITS
    # =========================================================================

    if (!"unit" %in% class(page$width)) {
        if (!is.numeric(page$width)) {
            stop("`width` is neither a unit object nor a numeric value. ",
                "Cannot create `plotgardener` page.", call. = FALSE)
        }

        if (is.null(page$default.units)) {
            stop("`width` detected as numeric. `default.units` ",
                "must be specified.", call. = FALSE)
        }

        if (!page$default.units %in% validUnits) {
            message <- paste(c("Invalid default units. ",
                            "Options for page units are:", validUnits),
                collapse = " "
            )
            stop(message, call. = FALSE)
        }

        page$width <- unit(page$width, page$default.units)
    }

    if (!"unit" %in% class(page$height)) {
        if (!is.numeric(page$height)) {
            stop("`height` is neither a unit object or a numeric value. ",
                "Cannot create `plotgardener` page.", call. = FALSE)
        }

        if (is.null(page$default.units)) {
            stop("`height` detected as numeric. `default.units` ",
                "must be specified.",
                call. = FALSE
            )
        }


        if (!page$default.units %in% validUnits) {
            message <- paste(c("Invalid default units. Options for page ",
                            "units are:", validUnits),
                collapse = " "
            )
            stop(message, call. = FALSE)
        }

        page$height <- unit(page$height, page$default.units)
    }

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    errorcheck_pageCreate(object = page)

    # =========================================================================
    # ASSIGN PAGE PARAMETERS
    # =========================================================================
    page_height <- as.numeric(page$height)
    assign("page_height", page_height, envir = pgEnv)

    page_units <- gsub("[0-9]|[.]", "", page$width)
    assign("page_units", page_units, envir = pgEnv)

    page_width <- as.numeric(page$width)

    # =========================================================================
    # VIEWPORT
    # =========================================================================

    page_vp <- viewport(
        x = unit(0.5, "npc"), y = unit(0.5, "npc"),
        width = page$width, height = page$height,
        xscale = c(0, page_width),
        yscale = rev(c(0, page_height)),
        name = "page"
    )

    # =========================================================================
    # INITIALIZE GTREE
    # =========================================================================

    assign("guide_grobs", gTree(name = "guide_grobs", vp = page_vp),
        envir = pgEnv
    )
    
    # =========================================================================
    # BACKGROUND
    # =========================================================================
    
    background <- rectGrob(gp = gpar(fill = page$bg, col = NA))
    assign("guide_grobs",
            addGrob(
                gTree = get("guide_grobs", envir = pgEnv),
                child = background
            ),
            envir = pgEnv
    )
    # =========================================================================
    # SHOW OUTLINE/RULER/UNITS
    # =========================================================================

    if (page$showGuides == TRUE) {
        border <- rectGrob()
        assign("guide_grobs",
            addGrob(
                gTree = get("guide_grobs", envir = pgEnv),
                child = border
            ),
            envir = pgEnv
        )


        if (page_units == "inches") {
            div <- 1 / 16
            x0 <- 0
            x1 <- -1 / 32
            y0 <- page_height
            y1 <- page_height + 1 / 32
            xsegs <- segmentsGrob(
                x0 = seq(0, page_width, div), y0 = y0,
                x1 = seq(0, page_width, div), y1 = y1,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            ysegs <- segmentsGrob(
                x0 = x0, y0 = seq(0, page_height, div),
                x1 = x1, y1 = seq(0, page_height, div),
                default.units = "native",
                gp = gpar(col = "black")
            )

            for (i in seq(1, 4)) {
                div <- div * 2
                x0 <- x1
                x1 <- x1 - 1 / 32
                y0 <- y1
                y1 <- y1 + 1 / 32
                v <- segmentsGrob(
                    x0 = xsegs$x0[xsegs$x0 %in% unit(seq(0, page_width, div), 
                                                     page_units)],
                    y0 = y0,
                    x1 = xsegs$x0[xsegs$x0 %in% unit(seq(0, page_width, div),
                                                     page_units)],
                    y1 = y1, default.units = page_units,
                    gp = gpar(col = "black")
                )
                h <- segmentsGrob(
                    x0 = x0,
                    y0 = ysegs$y1[ysegs$y1 %in% unit(seq(0, page_height, div),
                                                     "native")],
                    x1 = x1,
                    y1 = ysegs$y1[ysegs$y1 %in% unit(seq(0, page_height, div),
                                                     "native")],
                    default.units = "native",
                    gp = gpar(col = "black")
                )

                assign("guide_grobs",
                    addGrob(
                        gTree = get("guide_grobs", envir = pgEnv),
                        child = v
                    ),
                    envir = pgEnv
                )
                assign("guide_grobs",
                    addGrob(
                        gTree = get("guide_grobs", envir = pgEnv),
                        child = h
                    ),
                    envir = pgEnv
                )
            }

            hLabel <- textGrob(
                label = seq(0, page_width, div),
                x = seq(0, page_width, div),
                y = y0, vjust = -0.5,
                default.units = page_units
            )
            vLabel <- textGrob(
                label = seq(0, page_height, div),
                x = x0,
                y = seq(0, page_height, div),
                hjust = 1.5,
                default.units = "native"
            )
        } else if (page_units == "cm") {
            tickH <- convertUnit(
                x = unit(1 / 32, "inches"),
                unitTo = "cm", valueOnly = TRUE
            )
            div <- 1 / 10
            x0 <- 0
            x1 <- -3 * tickH
            y0 <- page_height
            y1 <- page_height + 3 * tickH
            xsegs <- segmentsGrob(
                x0 = seq(0, page_width, div),
                y0 = y0,
                x1 = seq(0, page_width, div),
                y1 = y1,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            ysegs <- segmentsGrob(
                x0 = x0,
                y0 = seq(0, page_height, div),
                x1 = x1,
                y1 = seq(0, page_height, div),
                default.units = "native",
                gp = gpar(col = "black")
            )


            div2 <- 1 / 2
            v2 <- segmentsGrob(
                x0 = seq(0, page_width, div2),
                y0 = page_height + 3 * tickH,
                x1 = seq(0, page_width, div2),
                y1 = page_height + 4 * tickH,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            h2 <- segmentsGrob(
                x0 = -3 * tickH,
                y0 = seq(0, page_height, div2),
                x1 = -4 * tickH,
                y1 = seq(0, page_height, div2),
                default.units = "native",
                gp = gpar(col = "black")
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = v2
                ),
                envir = pgEnv
            )
            assign("guide_grobs",
                addGrob(gTree = get("guide_grobs",
                    envir = pgEnv
                ), child = h2),
                envir = pgEnv
            )

            div3 <- 1
            v3 <- segmentsGrob(
                x0 = seq(0, page_width, div3),
                y0 = page_height + 4 * tickH,
                x1 = seq(0, page_width, div3),
                y1 = page_height + 5 * tickH,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            h3 <- segmentsGrob(
                x0 = -4 * tickH,
                y0 = seq(0, page_height, div3),
                x1 = -5 * tickH,
                y1 = seq(0, page_height, div3),
                default.units = "native",
                gp = gpar(col = "black")
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = v3
                ),
                envir = pgEnv
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = h3
                ),
                envir = pgEnv
            )

            hLabel <- textGrob(
                label = seq(0, page_width, div3),
                x = seq(0, page_width, div3),
                y = page_height + 4 * tickH,
                vjust = -0.5,
                default.units = page_units
            )
            vLabel <- textGrob(
                label = seq(0, page_height, div3),
                x = -4 * tickH,
                y = seq(0, page_height, div3),
                hjust = 1.5,
                default.units = "native"
            )
        } else if (page_units == "mm") {
            tickH <- convertUnit(
                x = unit(1 / 32, "inches"),
                unitTo = "mm", valueOnly = TRUE
            )
            div <- 1
            x0 <- 0
            x1 <- -3 * tickH
            y0 <- page_height
            y1 <- page_height + 3 * tickH
            xsegs <- segmentsGrob(
                x0 = seq(0, page_width, div),
                y0 = y0,
                x1 = seq(0, page_width, div),
                y1 = y1,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            ysegs <- segmentsGrob(
                x0 = x0,
                y0 = seq(0, page_height, div),
                x1 = x1,
                y1 = seq(0, page_height, div),
                default.units = "native",
                gp = gpar(col = "black")
            )

            div2 <- 5
            v2 <- segmentsGrob(
                x0 = seq(0, page_width, div2),
                y0 = page_height + 3 * tickH,
                x1 = seq(0, page_width, div2),
                y1 = page_height + 4 * tickH,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            h2 <- segmentsGrob(
                x0 = -3 * tickH,
                y0 = seq(0, page_height, div2),
                x1 = -4 * tickH,
                y1 = seq(0, page_height, div2),
                default.units = "native",
                gp = gpar(col = "black")
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = v2
                ),
                envir = pgEnv
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = h2
                ),
                envir = pgEnv
            )

            div3 <- 10
            v3 <- segmentsGrob(
                x0 = seq(0, page_width, div3),
                y0 = page_height + 4 * tickH,
                x1 = seq(0, page_width, div3),
                y1 = page_height + 5 * tickH,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            h3 <- segmentsGrob(
                x0 = -4 * tickH,
                y0 = seq(0, page_height, div3),
                x1 = -5 * tickH,
                y1 = seq(0, page_height, div3),
                default.units = "native",
                gp = gpar(col = "black")
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = v3
                ),
                envir = pgEnv
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = h3
                ),
                envir = pgEnv
            )

            hLabel <- textGrob(
                label = seq(0, page_width, div3),
                x = seq(0, page_width, div3),
                y = page_height + 4 * tickH,
                vjust = -0.5, default.units = page_units
            )
            vLabel <- textGrob(
                label = seq(0, page_height, div3),
                x = -4 * tickH,
                y = seq(0, page_height, div3),
                hjust = 1.5, default.units = "native"
            )
        } else {
            tickH <- convertUnit(
                x = unit(1 / 32, "inches"),
                unitTo = page_units, valueOnly = TRUE
            )
            div <- 1 / 10
            x0 <- 0
            x1 <- -3 * tickH
            y0 <- page_height
            y1 <- page_height + 3 * tickH
            xsegs <- segmentsGrob(
                x0 = seq(0, page_width, div),
                y0 = y0,
                x1 = seq(0, page_width, div),
                y1 = y1,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            ysegs <- segmentsGrob(
                x0 = x0,
                y0 = seq(0, page_height, div),
                x1 = x1,
                y1 = seq(0, page_height, div),
                default.units = "native",
                gp = gpar(col = "black")
            )


            div2 <- 1 / 2
            v2 <- segmentsGrob(
                x0 = seq(0, page_width, div2),
                y0 = page_height + 3 * tickH,
                x1 = seq(0, page_width, div2),
                y1 = page_height + 4 * tickH,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            h2 <- segmentsGrob(
                x0 = -3 * tickH,
                y0 = seq(0, page_height, div2),
                x1 = -4 * tickH,
                y1 = seq(0, page_height, div2),
                default.units = "native",
                gp = gpar(col = "black")
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = v2
                ),
                envir = pgEnv
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = h2
                ),
                envir = pgEnv
            )

            div3 <- 1
            v3 <- segmentsGrob(
                x0 = seq(0, page_width, div3),
                y0 = page_height + 4 * tickH,
                x1 = seq(0, page_width, div3),
                y1 = page_height + 5 * tickH,
                default.units = page_units,
                gp = gpar(col = "black")
            )
            h3 <- segmentsGrob(
                x0 = -4 * tickH,
                y0 = seq(0, page_height, div3),
                x1 = -5 * tickH,
                y1 = seq(0, page_height, div3),
                default.units = "native",
                gp = gpar(col = "black")
            )
            assign("guide_grobs",
                addGrob(
                    gTree = get("guide_grobs", envir = pgEnv),
                    child = v3
                ),
                envir = pgEnv
            )
            assign("guide_grobs",
                addGrob(gTree = get("guide_grobs",
                    envir = pgEnv
                ), child = h3),
                envir = pgEnv
            )

            hLabel <- textGrob(
                label = seq(0, page_width, div3),
                x = seq(0, page_width, div3),
                y = page_height + 4 * tickH,
                vjust = -0.5, default.units = page_units
            )
            vLabel <- textGrob(
                label = seq(0, page_height, div3),
                x = -4 * tickH,
                y = seq(0, page_height, div3),
                hjust = 1.5, default.units = "native"
            )
        }

        assign("guide_grobs",
            addGrob(
                gTree = get("guide_grobs", envir = pgEnv),
                child = xsegs
            ),
            envir = pgEnv
        )
        assign("guide_grobs",
            addGrob(
                gTree = get("guide_grobs", envir = pgEnv),
                child = ysegs
            ),
            envir = pgEnv
        )
        assign("guide_grobs",
            addGrob(
                gTree = get("guide_grobs", envir = pgEnv),
                child = hLabel
            ),
            envir = pgEnv
        )
        assign("guide_grobs",
            addGrob(
                gTree = get("guide_grobs", envir = pgEnv),
                child = vLabel
            ),
            envir = pgEnv
        )


        ## Unit annotation
        unitLabel <- textGrob(
            label = page_units, x = 0,
            y = page_height, hjust = 1.75,
            vjust = -1.5, default.units = page_units,
            just = c("right", "bottom")
        )
        assign("guide_grobs",
            addGrob(
                gTree = get("guide_grobs", envir = pgEnv),
                child = unitLabel
            ),
            envir = pgEnv
        )
        
        # =====================================================================
        # X AND Y GRIDLINES
        # =====================================================================
        
        ## Draw x and y gridlines
        tryCatch(expr = {
            xGrid <- segmentsGrob(
                x0 = seq(0, page_width, page$xgrid),
                y0 = 0,
                x1 = seq(0, page_width, page$xgrid),
                y1 = page$height,
                default.units = page_units,
                gp = gpar(col = "grey50", lty = 2, lwd = 0.5)
            )
            
            yGrid <- segmentsGrob(
                x0 = 0,
                y0 = seq(0, page_height, page$ygrid),
                x1 = page_width,
                y1 = seq(0, page_height, page$ygrid),
                default.units = "native",
                gp = gpar(col = "grey50", lty = 2, lwd = 0.5)
            )
            
            assign("guide_grobs",
                    addGrob(
                        gTree = get("guide_grobs", envir = pgEnv),
                        child = xGrid
                    ),
                    envir = pgEnv
            )
            assign("guide_grobs",
                    addGrob(
                        gTree = get("guide_grobs", envir = pgEnv),
                        child = yGrid
                    ),
                    envir = pgEnv
            )
        }, error = function(e) {
            return()
        })
        
        }

    # =========================================================================
    # DRAW GROBS
    # =========================================================================

    grid.draw(get("guide_grobs", envir = pgEnv))
    downViewport("page")
}
