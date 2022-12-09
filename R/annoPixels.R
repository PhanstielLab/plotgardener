#' Annotate pixels in a Hi-C plot
#' 
#' @usage annoPixels(
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
#' @param plot Hi-C plot object from \code{plotHicSquare} or
#' \code{plotHicTriangle} on which to annotate pixels.
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
#' \item{\code{"bottom"}: }{Pixels will be annotated on the bottom diagonal
#' half of a square Hi-C plot.}
#' }
#' @param shift Numeric specifying the number of pixels on either end of
#' main pixel in a box or circle. Numeric specifying number of pixels
#' for the length of an arrow.
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#' @param quiet A logical indicating whether or not to print messages.
#'
#' @return Returns a \code{pixel} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data and BEDPE data
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
#' data("IMR90_DNAloops_pairs")
#'
#' ## Create page
#' pageCreate(width = 4.5, height = 4, default.units = "inches")
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- plotHicSquare(
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
#' pixels <- annoPixels(
#'     plot = hicPlot, data = IMR90_DNAloops_pairs, type = "box",
#'     half = "both"
#' )
#'
#' ## Annotate loops on one side of Hi-C plot with arrows
#' ## and the other side with circles
#' pagePlotRemove(plot = pixels)
#' pixels1 <- annoPixels(
#'     plot = hicPlot, data = IMR90_DNAloops_pairs,
#'     type = "arrow", half = "top", shift = 8
#' )
#' pixels2 <- annoPixels(
#'     plot = hicPlot, data = IMR90_DNAloops_pairs,
#'     type = "circle", half = "bottom"
#' )
#'
#' ## Annotate heatmap legend
#' annoHeatmapLegend(
#'     plot = hicPlot,
#'     x = 3.6, y = 0.5, width = 0.12, height = 1.2,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' annoGenomeLabel(
#'     plot = hicPlot, x = 0.5, y = 3.53, scale = "Mb",
#'     just = c("left", "top")
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @export
annoPixels <- function(plot, data, type = "box", half = "inherit",
                        shift = 4, params = NULL, quiet = FALSE, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================
    ## Define a function to catch errors for annoPixels
    errorcheck_annoLoops <- function(hic, loops, half, type, quiet) {

        ###### hic #####

        ## check type of input for hic
        if (!class(hic) %in% c(
            "hicSquare", "hicTriangle",
            "hicRectangle"
        )) {
            stop("Input plot must be a plot of class \'hicSquare\', ",
                "\'hicTriangle\', or \'hicRectangle\'.", call. = FALSE)
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
        if (is(hic, "hicSquare")) {
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
        } else if (is(hic, "hicTriangle") |
            is(hic, "hicRectangle")) {
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
    subset_loops <- function(hic, loopData, object, plotObject) {

        ## chrom always in col1
        ## altchrom always in col4
        ## triangle hic plots will not have altchrom parameters
        if (is(hic, "hicTriangle")) {
            loops_subset <- loopData[which(
                loopData[, "chrom1"] == object$chrom &
                loopData[, "chrom2"] == object$chrom &
                loopData[, "start1"] >= object$chromstart &
                loopData[, "end1"] <= object$chromend &
                loopData[, "start2"] >= object$chromstart &
                loopData[, "end2"] <= object$chromend), ]
        } else if (is(hic, "hicRectangle")){
            loops_subset <- loopData[which(
                loopData[, "chrom1"] == object$chrom &
                    loopData[, "chrom2"] == object$chrom &
                    loopData[, "start1"] >= plotObject$chromstartAdjusted &
                    loopData[, "end1"] <= plotObject$chromendAdjusted &
                    loopData[, "start2"] >= plotObject$chromstartAdjusted &
                    loopData[, "end2"] <= plotObject$chromendAdjusted), ]
        } else {
            loops_subset <- loopData[which(
                loopData[, "chrom1"] == object$chrom &
                loopData[, "chrom2"] == object$altchrom &
                loopData[, "start1"] >= object$chromstart &
                loopData[, "end1"] <= object$chromend &
                loopData[, "start2"] >= object$altchromstart &
                loopData[, "end2"] <= object$altchromend), ]
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = rect1
                ),
                envir = pgEnv
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = rect1
                ),
                envir = pgEnv
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = rect1
                ),
                envir = pgEnv
            )
            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = rect2
                ),
                envir = pgEnv
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = circ1
                ),
                envir = pgEnv
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = circ1
                ),
                envir = pgEnv
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = circ1
                ),
                envir = pgEnv
            )
            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = circ2
                ),
                envir = pgEnv
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = arrow1
                ),
                envir = pgEnv
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = arrow1
                ),
                envir = pgEnv
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
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = arrow1
                ),
                envir = pgEnv
            )
            assign("loop_grobs",
                addGrob(
                    gTree = get("loop_grobs", envir = pgEnv),
                    child = arrow2
                ),
                envir = pgEnv
            )
        }
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================
    
    loopsInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "loopsInternal"
    )

    ## Set gp
    loopsInternal$gp <- setGP(
        gpList = gpar(),
        params = loopsInternal, ...
    )

    # =========================================================================
    # INITIALIZE OBJECT: GET REGION/DIMENSIONS FROM HIC PLOT INPUT
    # =========================================================================

    loops <- structure(list(
        chrom = loopsInternal$plot$chrom,
        chromstart = loopsInternal$plot$chromstart,
        chromend = loopsInternal$plot$chromend,
        altchrom = loopsInternal$plot$altchrom,
        altchromstart = loopsInternal$plot$altchromstart,
        altchromend = loopsInternal$plot$altchromend,
        assembly = loopsInternal$plot$assembly,
        x = loopsInternal$plot$x,
        y = loopsInternal$plot$y,
        width = loopsInternal$plot$width,
        height = loopsInternal$plot$height,
        just = loopsInternal$plot$just, grobs = NULL
    ),
    class = "pixel"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_page(error = "Cannot annotate Hi-C pixels without a
                `plotgardener` page.")
    if (is.null(loopsInternal$plot)) stop("argument \"plot\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(loopsInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)

    errorcheck_annoLoops(
        hic = loopsInternal$plot,
        loops = loopsInternal$data,
        half = loopsInternal$half,
        type = loopsInternal$type,
        quiet = loopsInternal$quiet
    )

    # =========================================================================
    # PARSE INHERITED HALF
    # =========================================================================

    half <- loopsInternal$half
    if (half == "inherit") {
        half <- inherit_half(hic = loopsInternal$plot)
    }

    if (is(loopsInternal$plot, "hicTriangle")  |
        is(loopsInternal$plot, "hicRectangle")) {
        if (loopsInternal$plot$flip == TRUE){
            half <- "bottom"
        } else {
            half <- "top"
        }
        
    }

    # =========================================================================
    # READ IN FILE, DATAFRAME OR GINTERACTIONS
    # =========================================================================

    loopData <- read_pairedData(data = loopsInternal$data,
                            assembly = loops$assembly,
                            warning = TRUE)

    ## chrom format and data chrom format
    chromDataAgreement(data = loopData, chrom = loops$chrom,
                    type = "pairs")

    # =========================================================================
    # SUBSET FOR LOOPS IN REGION
    # =========================================================================
    loops_subset <- subset_loops(
        hic = loopsInternal$plot, loopData = loopData,
        object = loops, plotObject = loopsInternal$plot
    )
    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "pixel",
        length(grep(
            pattern = "pixel",
            x = currentViewports
        )) + 1
    )

    ## Make viewport based on hic input viewport
    if (is(loopsInternal$plot, "hicSquare")) {
        vp <- viewport(
            height = loopsInternal$plot$grobs$vp$height,
            width = loopsInternal$plot$grobs$vp$width,
            x = loopsInternal$plot$grobs$vp$x,
            y = loopsInternal$plot$grobs$vp$y,
            clip = "on",
            xscale = loopsInternal$plot$grobs$vp$xscale,
            yscale = loopsInternal$plot$grobs$vp$yscale,
            just = loopsInternal$plot$grobs$vp$justification,
            name = vp_name
        )
    } else if (is(loopsInternal$plot, "hicTriangle")) {
        vp <- viewport(
            height = loopsInternal$plot$grobs$vp$height,
            width = loopsInternal$plot$grobs$vp$width,
            x = unit(0, "npc"),
            y =  unit(0, "npc"),
            xscale = loopsInternal$plot$grobs$vp$xscale,
            yscale = loopsInternal$plot$grobs$vp$yscale,
            just = c("left", "bottom"),
            name = vp_name,
            angle = -45
        )
        
        if (loopsInternal$plot$flip == TRUE){
            vp$y <- unit(1, "npc")

        }
        
    } else if (is(loopsInternal$plot, "hicRectangle")) {
        side <- convertUnit(loopsInternal$plot$grobs$vp$width,
            unitTo = get("page_units", pgEnv)
        )

        vp <- viewport(
            height = side, width = side,
            x = unit(loopsInternal$plot$grobs$vp$xscale[1], "native"), 
            y = unit(0, "npc"),
            xscale = loopsInternal$plot$grobs$vp$xscale,
            yscale = loopsInternal$plot$grobs$vp$yscale,
            just = c("left", "bottom"),
            name = vp_name,
            angle = -45
        )
        if (loopsInternal$plot$flip == TRUE){
            vp$y <- unit(1, "npc")
            
        }
    }

    # =========================================================================
    # INITIALIZE GTREE OF GROBS
    # =========================================================================

    assign("loop_grobs", gTree(vp = vp), envir = pgEnv)

    # =========================================================================
    # PLOT
    # =========================================================================

    if (nrow(loops_subset) > 0) {
        if (loopsInternal$type == "box") {
            loopsInternal$gp$fill <- NA
            invisible(apply(loops_subset, 1, boxAnnotation,
                hic = loopsInternal$plot, object = loopsInternal,
                shift = loopsInternal$shift, half = half
            ))
        } else if (loopsInternal$type == "circle") {
            loopsInternal$gp$fill <- NA
            invisible(apply(loops_subset, 1, circleAnnotation,
                hic = loopsInternal$plot, object = loopsInternal,
                shift = loopsInternal$shift, half = half
            ))
        } else if (loopsInternal$type == "arrow") {
            if (is.null(loopsInternal$gp$col) &
                is.null(loopsInternal$gp$fill)) {
                loopsInternal$gp$fill <- "black"
            } else {
                if (is.null(loopsInternal$gp$fill)) {
                    loopsInternal$gp$fill <- loopsInternal$gp$col
                }
            }

            invisible(apply(loops_subset, 1, arrowAnnotation,
                hic = loopsInternal$plot, object = loopsInternal,
                shift = loopsInternal$shift, half = half
            ))
        }
    } else {
        warning("No pixels found in region.", call. = FALSE)
    }


    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    loops$grobs <- get("loop_grobs", envir = pgEnv)
    
    if (is(loopsInternal$plot, "hicRectangle") | 
        is(loopsInternal$plot, "hicTriangle")){
        seekViewport(name = loopsInternal$plot$outsideVP$name)
        grid.draw(loops$grobs)
        seekViewport("page")
    } else {
        grid.draw(loops$grobs)
    }
    
    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("pixel[", vp_name, "]")
    invisible(loops)
}
