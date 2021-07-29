#' Annotate pixels in a Hi-C plot
#' 
#' @usage bbAnnoPixels(
#'     plot,
#'     data,
#'     type = "box",
#'     half = "inherit",
#'     shift = 4,
#'     params = NULL,
#'     quiet = FALSE,
#'     ...
#' )
#'
#' @param plot Hi-C plot object from \code{bbPlotHicSquare} or
#' \code{bbPlotHicTriangle} on which to annotate pixels.
#' @param data A string specifying the BEDPE file path, a dataframe in BEDPE
#' format specifying pixel positions, or a
#' \link[InteractionSet]{GInteractions} object specifying pixel
#' positions.
#' @param type Character value specifying type of annotation.
#' Default value is \code{type = "box"}. Options are:
#' \itemize{
#' \item{\code{"box"}: }{Boxes are drawn around each pixel.}
#' \item{\code{"circle"}: }{Circles are drawn around each pixel.}
#' \item{\code{"arrow"}: }{Arrows are drawn pointing to each pixel.}
#' }
#' @param half Character value specifying which half of hic plots
#' to annotate. Triangle Hi-C plots will always default to the entirety of
#' the triangular plot. Default value is \code{half = "inherit"}. Options are:
#' \itemize{
#' \item{\code{"inherit"}: }{Pixels will be annotated on the \code{half}
#' inherited by the input Hi-C plot.}
#' \item{\code{"both"}: }{Pixels will be annotated on both halves of the
#' diagonal of a square Hi-C plot.}
#' \item{\code{"top"}: }{Pixels will be annotated on the upper diagonal
#' half of a square Hi-C plot.}
#' \item{\code{"bottom"}: }{Pixels will be annotated ont the bottom diagonal
#' half of a square Hi-C plot.}
#' }
#' @param shift Numeric specifying the number of pixels on either end of
#' main pixel in a box or circle. Numeric specifying number of pixels
#' for the length of an arrow.
#' @param params An optional \link[BentoBox]{bbParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#' @param quiet A logical indicating whether or not to print messages.
#'
#' @return Returns a \code{bb_pixel} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data and BEDPE data
#' library(BentoBoxData)
#' data("IMR90_HiC_10kb")
#' data("IMR90_DNAloops_pairs")
#'
#' ## Create BentoBox page
#' bbPageCreate(width = 4.5, height = 4, default.units = "inches")
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- bbPlotHicSquare(
#'     data = IMR90_HiC_10kb, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     x = 0.5, y = 0.5, width = 3, height = 3,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate loops of both sides of Hi-C plot with squares
#' pixels <- bbAnnoPixels(
#'     plot = hicPlot, data = IMR90_DNAloops_pairs, type = "box",
#'     half = "both"
#' )
#'
#' ## Annotate loops on one side of Hi-C plot with arrows
#' ## and the other side with circles
#' bbPagePlotRemove(plot = pixels)
#' pixels1 <- bbAnnoPixels(
#'     plot = hicPlot, data = IMR90_DNAloops_pairs,
#'     type = "arrow", half = "top", shift = 8
#' )
#' pixels2 <- bbAnnoPixels(
#'     plot = hicPlot, data = IMR90_DNAloops_pairs,
#'     type = "circle", half = "bottom"
#' )
#'
#' ## Annotate heatmap legend
#' bbAnnoHeatmapLegend(
#'     plot = hicPlot,
#'     x = 3.6, y = 0.5, width = 0.12, height = 1.2,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bbAnnoGenomeLabel(
#'     plot = hicPlot, x = 0.5, y = 3.53, scale = "Mb",
#'     just = c("left", "top")
#' )
#'
#' ## Hide page guides
#' bbPageGuideHide()
#' @export
bbAnnoPixels <- function(plot, data, type = "box", half = "inherit",
                        shift = 4, params = NULL, quiet = FALSE, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================
    ## Define a function to catch errors for bbAnnoPixels
    errorcheck_bbAnnoLoops <- function(hic, loops, half, type, quiet) {

        ###### hic #####

        ## check type of input for hic
        if (!class(hic) %in% c(
            "bb_hicSquare", "bb_hicTriangle",
            "bb_hicRectangle"
        )) {
            stop("Input plot must be a plot of class \'bb_hicSquare\', ",
                "\'bb_hicTriangle\', or \'bb_hicRectangle\'.", call. = FALSE)
        }

        ###### loops #####

        ## if data.frame/data.table needs to be properly formatted
        if ("data.frame" %in% class(loops) && ncol(loops) < 6) {
            stop("Invalid dataframe format. ",
                "Dataframe must be in BEDPE format.", call. = FALSE)
        }

        if ("data.frame" %in% class(loops) && nrow(loops) < 1) {
            stop("\'data\' input contains no values.", call. = FALSE)
        }


        ## if it's a file path, it needs to exist
        if (!"data.frame" %in% class(loops)) {
            if (!is(loops, "GInteractions")) {
                ## File existence
                if (!file.exists(loops)) {
                    stop("File", loops, "does not exist.", call. = FALSE)
                }
            }
        }

        ###### half #####

        ## half needs to be a valid option
        if (!half %in% c("inherit", "both", "top", "bottom")) {
            stop("Invalid \'half\'.  Options are \'inherit\',
                \'both\', \'top\', or \'bottom\'.", call. = FALSE)
        }

        ## half needs to be able to align with what kind of hic plot is plotted
        if (is(hic, "bb_hicSquare")) {
            if (hic$chrom == hic$altchrom) {
                if ((hic$half == "top" | hic$half == "bottom") &&
                    (half == "both")) {
                    stop("Invalid \'half\' of plot to annotate.",
                        call. = FALSE
                    )
                }

                if (hic$half == "top" & half == "bottom") {
                    stop("Invalid \'half\' of plot to annotate.",
                        call. = FALSE
                    )
                }

                if (hic$half == "bottom" & half == "top") {
                    stop("Invalid \'half\' of plot to annotate.",
                        call. = FALSE
                    )
                }
            } else {
                if (hic$half == "bottom") {
                    if (!quiet) {
                        message("Attempting to annotate pixels where",
                            hic$chrom, "is on the x-axis and",
                            hic$altchrom, "is on the y-axis.",
                            call. = FALSE
                        )
                    }
                } else if (hic$half == "top") {
                    if (!quiet) {
                        message("Attempting to annotate pixels where",
                            hic$altchrom, "is on the x-axis and",
                            hic$chrom, "is on the y-axis.",
                            call. = FALSE
                        )
                    }
                }
            }
        } else if (is(hic, "bb_hicTriangle") |
            is(hic, "bb_hicRectangle")) {
            if (half == "both" | half == "bottom") {
                warning("Plot of class \'",
                    class(hic),
                    "\' detected. Pixels will automatically be annotated ",
                    "in the upper triangular of the plot.",
                    call. = FALSE
                )
            }
        }

        ###### annotation #####

        ## Check type of annotation
        if (!type %in% c("box", "circle", "arrow")) {
            stop("Invalid \'type\' of annotation.  Options are \'box\', ",
                "\'circle\', or \'arrow\'.", call. = FALSE)
        }
    }

    ## Define a function that subsets loop data for hic region
    subset_loops <- function(hic, loops, object) {

        ## chrom always in col1
        ## altchrom always in col4
        ## triangle hic plots will not have altchrom parameters
        if (is(hic, "bb_hicTriangle") | is(hic, "bb_hicRectangle")) {
            loops_subset <- loops[which(loops[, "chrom1"] == object$chrom &
                loops[, "chrom2"] == object$chrom &
                loops[, "start1"] >= object$chromstart &
                loops[, "end1"] <= object$chromend &
                loops[, "start2"] >= object$chromstart &
                loops[, "end2"] <= object$chromend), ]
        } else {
            loops_subset <- loops[which(loops[, "chrom1"] == object$chrom &
                loops[, "chrom2"] == object$altchrom &
                loops[, "start1"] >= object$chromstart &
                loops[, "end1"] <= object$chromend &
                loops[, "start2"] >= object$altchromstart &
                loops[, "end2"] <= object$altchromend), ]
        }

        return(loops_subset)
    }

    ## Define a function to add box annotation
    boxAnnotation <- function(df, hic, object, shift, half) {
        side <- (utils::type.convert(df["end2"], as.is = TRUE) - 
                    utils::type.convert(df["start2"], as.is = TRUE)) +
            (2 * shift * hic$resolution)

        if (half == "bottom") {
            center_x <- 0.5 * (utils::type.convert(df["start2"], 
                                        as.is = TRUE) + 
                        utils::type.convert(df["end2"], 
                                        as.is = TRUE))
            center_y <- 0.5 * (utils::type.convert(df["start1"], 
                                        as.is = TRUE) + 
                        utils::type.convert(df["end1"], 
                                        as.is = TRUE))
            rect1 <- rectGrob(
                x = center_x, y = center_y, width = side,
                height = side, default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = rect1
                ),
                envir = bbEnv
            )
        } else if (half == "top") {
            center_x <- 0.5 * (utils::type.convert(df["start1"], as.is = TRUE) 
                            + utils::type.convert(df["end1"], as.is = TRUE))
            center_y <- 0.5 * (utils::type.convert(df["start2"], as.is = TRUE) 
                            + utils::type.convert(df["end2"], as.is = TRUE))
            rect1 <- rectGrob(
                x = center_x, y = center_y, width = side,
                height = side, default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = rect1
                ),
                envir = bbEnv
            )
        } else if (half == "both") {

            ## BOTTOM
            center_x1 <- 0.5 * (utils::type.convert(df["start2"], as.is = TRUE) 
                                + utils::type.convert(df["end2"], as.is = TRUE))
            center_y1 <- 0.5 * (utils::type.convert(df["start1"], as.is = TRUE) 
                                + utils::type.convert(df["end1"], as.is = TRUE))

            ## TOP
            center_x2 <- 0.5 * (utils::type.convert(df["start1"], as.is = TRUE) 
                                + utils::type.convert(df["end1"], as.is = TRUE))
            center_y2 <- 0.5 * (utils::type.convert(df["start2"], as.is = TRUE) 
                                + utils::type.convert(df["end2"], as.is = TRUE))

            rect1 <- rectGrob(
                x = center_x1, y = center_y1, width = side,
                height = side, default.units = "native",
                gp = object$gp
            )
            rect2 <- rectGrob(
                x = center_x2, y = center_y2, width = side,
                height = side, default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = rect1
                ),
                envir = bbEnv
            )
            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = rect2
                ),
                envir = bbEnv
            )
        }
    }

    ## Define a function to add circle annotation
    circleAnnotation <- function(df, hic, object, shift, half) {
        radius <- (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) - 
                        utils::type.convert(df["start2"], 
                                as.is = TRUE))) +
            (shift * hic$resolution)

        if (half == "bottom") {
            center_x <- 0.5 * (utils::type.convert(df["start2"], as.is = TRUE)
                            + utils::type.convert(df["end2"], as.is = TRUE))
            center_y <- 0.5 * (utils::type.convert(df["start1"], as.is = TRUE) 
                            + utils::type.convert(df["end1"], as.is = TRUE))
            circ1 <- circleGrob(
                x = center_x, y = center_y,
                r = radius, default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = circ1
                ),
                envir = bbEnv
            )
        } else if (half == "top") {
            center_x <- 0.5 * (utils::type.convert(df["start1"], as.is = TRUE) 
                            + utils::type.convert(df["end1"], as.is = TRUE))
            center_y <- 0.5 * (utils::type.convert(df["start2"], as.is = TRUE) 
                            + utils::type.convert(df["end2"], as.is = TRUE))
            circ1 <- circleGrob(
                x = center_x, y = center_y,
                r = radius, default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = circ1
                ),
                envir = bbEnv
            )
        } else if (half == "both") {

            ## BOTTOM
            center_x1 <- 0.5 * (utils::type.convert(df["start2"], 
                                                    as.is = TRUE) + 
                            utils::type.convert(df["end2"], 
                                            as.is = TRUE))
            center_y1 <- 0.5 * (utils::type.convert(df["start1"], 
                                                    as.is = TRUE) + 
                            utils::type.convert(df["end1"], 
                                            as.is = TRUE))

            ## TOP
            center_x2 <- 0.5 * (utils::type.convert(df["start1"], as.is = TRUE) 
                                + utils::type.convert(df["end1"], as.is = TRUE))
            center_y2 <- 0.5 * (utils::type.convert(df["start2"], as.is = TRUE) 
                                + utils::type.convert(df["end2"], as.is = TRUE))

            circ1 <- circleGrob(
                x = center_x1, y = center_y1,
                r = radius, default.units = "native",
                gp = object$gp
            )
            circ2 <- circleGrob(
                x = center_x2, y = center_y2,
                r = radius, default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = circ1
                ),
                envir = bbEnv
            )
            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = circ2
                ),
                envir = bbEnv
            )
        }
    }

    ## Define a function to add arrow annotation
    arrowAnnotation <- function(df, hic, object, shift, half) {
        if (half == "bottom") {
            x0 <- utils::type.convert(df["end2"], as.is = TRUE) + 
                (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) -
                utils::type.convert(df["start2"], as.is = TRUE)))
            y0 <- utils::type.convert(df["start1"], as.is = TRUE) - 
                (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) -
                utils::type.convert(df["start2"], as.is = TRUE)))

            arrow1 <- segmentsGrob(
                x0 = x0, y0 = y0,
                x1 = x0 + (shift * hic$resolution),
                y1 = y0 - (shift * hic$resolution),
                arrow = arrow(
                    length = unit(0.1, "inches"),
                    ends = "first",
                    type = "closed"
                ),
                default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = arrow1
                ),
                envir = bbEnv
            )
        } else if (half == "top") {
            x0 <- utils::type.convert(df["start1"], as.is = TRUE) - 
                (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) -
                utils::type.convert(df["start2"], as.is = TRUE)))
            y0 <- utils::type.convert(df["end2"], as.is = TRUE) +
                (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) -
                utils::type.convert(df["start2"], as.is = TRUE)))

            arrow1 <- segmentsGrob(
                x0 = x0, y0 = y0,
                x1 = x0 - (shift * hic$resolution),
                y1 = y0 + (shift * hic$resolution),
                arrow = arrow(
                    length = unit(0.1, "inches"),
                    ends = "first",
                    type = "closed"
                ),
                default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = arrow1
                ),
                envir = bbEnv
            )
        } else if (half == "both") {

            ## BOTTOM
            x01 <- utils::type.convert(df["end2"], as.is = TRUE) + 
                (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) -
                utils::type.convert(df["start2"], as.is = TRUE)))
            y01 <- utils::type.convert(df["start1"], as.is = TRUE) - 
                (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) -
                utils::type.convert(df["start2"], as.is = TRUE)))

            ## TOP
            x02 <- utils::type.convert(df["start1"], as.is = TRUE) - 
                (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) -
                utils::type.convert(df["start2"], as.is = TRUE)))
            y02 <- utils::type.convert(df["end2"], as.is = TRUE) +
                (0.5 * (utils::type.convert(df["end2"], as.is = TRUE) -
                utils::type.convert(df["start2"], as.is = TRUE)))

            arrow1 <- segmentsGrob(
                x0 = x01, y0 = y01,
                x1 = x01 + (shift * hic$resolution),
                y1 = y01 - (shift * hic$resolution),
                arrow = arrow(
                    length = unit(0.1, "inches"),
                    ends = "first",
                    type = "closed"
                ),
                default.units = "native",
                gp = object$gp
            )
            arrow2 <- segmentsGrob(
                x0 = x02, y0 = y02,
                x1 = x02 - (shift * hic$resolution),
                y1 = y02 + (shift * hic$resolution),
                arrow = arrow(
                    length = unit(0.1, "inches"),
                    ends = "first",
                    type = "closed"
                ),
                default.units = "native",
                gp = object$gp
            )

            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = arrow1
                ),
                envir = bbEnv
            )
            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = bbEnv),
                    child = arrow2
                ),
                envir = bbEnv
            )
        }
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================
    
    bb_loopsInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_loopsInternal"
    )

    ## Set gp
    bb_loopsInternal$gp <- setGP(
        gpList = gpar(),
        params = bb_loopsInternal, ...
    )

    # =========================================================================
    # INITIALIZE OBJECT: GET REGION/DIMENSIONS FROM HIC PLOT INPUT
    # =========================================================================

    bb_loops <- structure(list(
        chrom = bb_loopsInternal$plot$chrom,
        chromstart = bb_loopsInternal$plot$chromstart,
        chromend = bb_loopsInternal$plot$chromend,
        altchrom = bb_loopsInternal$plot$altchrom,
        altchromstart = bb_loopsInternal$plot$altchromstart,
        altchromend = bb_loopsInternal$plot$altchromend,
        assembly = bb_loopsInternal$plot$assembly,
        x = bb_loopsInternal$plot$x,
        y = bb_loopsInternal$plot$y,
        width = bb_loopsInternal$plot$width,
        height = bb_loopsInternal$plot$height,
        just = bb_loopsInternal$plot$just, grobs = NULL
    ),
    class = "bb_pixel"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_bbpage(error = "Cannot annotate Hi-C pixels without a
                BentoBox page.")
    if (is.null(bb_loopsInternal$plot)) stop("argument \"plot\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(bb_loopsInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)

    errorcheck_bbAnnoLoops(
        hic = bb_loopsInternal$plot,
        loops = bb_loopsInternal$data,
        half = bb_loopsInternal$half,
        type = bb_loopsInternal$type,
        quiet = bb_loopsInternal$quiet
    )

    # =========================================================================
    # PARSE INHERITED HALF
    # =========================================================================

    half <- bb_loopsInternal$half
    if (half == "inherit") {
        half <- inherit_half(hic = bb_loopsInternal$plot)
    }

    if (is(bb_loopsInternal$plot, "bb_hicTriangle")  |
        is(bb_loopsInternal$plot, "bb_hicRectangle")) {
        half <- "top"
    }

    # =========================================================================
    # READ IN FILE, DATAFRAME OR GINTERACTIONS
    # =========================================================================

    loops <- read_pairedData(data = bb_loopsInternal$data,
                            assembly = bb_loops$assembly,
                            warning = TRUE)

    ## chrom format and data chrom format
    chromDataAgreement(data = loops, chrom = bb_loops$chrom,
                    type = "pairs")

    # =========================================================================
    # SUBSET FOR LOOPS IN REGION
    # =========================================================================
    loops_subset <- subset_loops(
        hic = bb_loopsInternal$plot, loops = loops,
        object = bb_loops
    )
    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_pixel",
        length(grep(
            pattern = "bb_pixel",
            x = currentViewports
        )) + 1
    )

    ## Make viewport based on hic input viewport
    if (is(bb_loopsInternal$plot, "bb_hicSquare")) {
        vp <- viewport(
            height = bb_loopsInternal$plot$grobs$vp$height,
            width = bb_loopsInternal$plot$grobs$vp$width,
            x = bb_loopsInternal$plot$grobs$vp$x,
            y = bb_loopsInternal$plot$grobs$vp$y,
            clip = "on",
            xscale = bb_loopsInternal$plot$grobs$vp$xscale,
            yscale = bb_loopsInternal$plot$grobs$vp$yscale,
            just = bb_loopsInternal$plot$grobs$vp$justification,
            name = vp_name
        )
    } else if (is(bb_loopsInternal$plot, "bb_hicTriangle")) {
        width <- convertUnit(bb_loopsInternal$plot$outsideVP$width,
            unitTo = get("page_units", bbEnv), valueOnly = TRUE
        )

        vp <- viewport(
            height = unit(width / sqrt(2), get("page_units", bbEnv)),
            width = unit(width / sqrt(2), get("page_units", bbEnv)),
            x = bb_loopsInternal$plot$outsideVP$x,
            y = bb_loopsInternal$plot$outsideVP$y,
            xscale = bb_loopsInternal$plot$grobs$vp$xscale,
            yscale = bb_loopsInternal$plot$grobs$vp$yscale,
            just = bb_loopsInternal$plot$outsideVP$justification,
            name = vp_name,
            angle = -45
        )
    } else if (is(bb_loopsInternal$plot, "bb_hicRectangle")) {
        side <- convertUnit(bb_loopsInternal$plot$grobs$vp$width,
            unitTo = get("page_units", bbEnv)
        )

        ## Get bottom left coord of outsideVP
        bottomLeft <- vp_bottomLeft(viewport = bb_loopsInternal$plot$outsideVP)

        ## Convert adjusted chromstart to page units within outsideVP
        ## and add to bottomLeft x
        seekViewport(name = bb_loopsInternal$plot$outsideVP$name)
        xCoord <- convertX(unit(bb_loopsInternal$plot$grobs$vp$x),
            unitTo = get("page_units", bbEnv)
        ) + bottomLeft[[1]]
        seekViewport(name = "bb_page")

        vp <- viewport(
            height = side, width = side,
            x = xCoord, y = bottomLeft[[2]],
            xscale = bb_loopsInternal$plot$grobs$vp$xscale,
            yscale = bb_loopsInternal$plot$grobs$vp$yscale,
            just = c("left", "bottom"),
            name = vp_name,
            angle = -45
        )
    }

    # =========================================================================
    # INITIALIZE GTREE OF GROBS
    # =========================================================================

    assign("loop_grobs", gTree(vp = vp), envir = bbEnv)

    # =========================================================================
    # PLOT
    # =========================================================================

    if (nrow(loops_subset) > 0) {
        if (bb_loopsInternal$type == "box") {
            bb_loopsInternal$gp$fill <- NA
            invisible(apply(loops_subset, 1, boxAnnotation,
                hic = bb_loopsInternal$plot, object = bb_loopsInternal,
                shift = bb_loopsInternal$shift, half = half
            ))
        } else if (bb_loopsInternal$type == "circle") {
            bb_loopsInternal$gp$fill <- NA
            invisible(apply(loops_subset, 1, circleAnnotation,
                hic = bb_loopsInternal$plot, object = bb_loopsInternal,
                shift = bb_loopsInternal$shift, half = half
            ))
        } else if (bb_loopsInternal$type == "arrow") {
            if (is.null(bb_loopsInternal$gp$col) &
                is.null(bb_loopsInternal$gp$fill)) {
                bb_loopsInternal$gp$fill <- "black"
            } else {
                if (is.null(bb_loopsInternal$gp$fill)) {
                    bb_loopsInternal$gp$fill <- bb_loopsInternal$gp$col
                }
            }

            invisible(apply(loops_subset, 1, arrowAnnotation,
                hic = bb_loopsInternal$plot, object = bb_loopsInternal,
                shift = bb_loopsInternal$shift, half = half
            ))
        }
    } else {
        warning("No pixels found in region.", call. = FALSE)
    }


    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    bb_loops$grobs <- get("loop_grobs", envir = bbEnv)
    grid.draw(bb_loops$grobs)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_pixel[", vp_name, "]")
    invisible(bb_loops)
}
