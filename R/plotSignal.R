#' Plot any kind of signal track data for a single chromosome
#' 
#' @usage plotSignal(
#'     data,
#'     binSize = NA,
#'     binCap = TRUE,
#'     negData = FALSE,
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
#'     linecolor = "#37a7db",
#'     fill = NA,
#'     ymax = 1,
#'     range = NULL,
#'     scale = FALSE,
#'     bg = NA,
#'     baseline = TRUE,
#'     baseline.color = "grey",
#'     baseline.lwd = 1,
#'     orientation = "h",
#'     x = NULL,
#'     y = NULL,
#'     width = NULL,
#'     height = NULL,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     draw = TRUE,
#'     params = NULL,
#'     ...
#' )
#'
#' @param data Data to be plotted as a character value specifying a
#' bigwig file path, a dataframe in BED format, or a
#' \link[GenomicRanges]{GRanges} object with metadata column \code{score}.
#' Either one \code{data} argument or a list of two can be provided, where
#' the second \code{data} will be plotted below the x-axis.
#' @param binSize A numeric specifying the length of each data
#' bin in basepairs. Default value is \code{binSize = NA}.
#' @param binCap A logical value indicating whether the function will
#' limit the number of data bins to 8,000.
#' Default value is \code{binCap = TRUE}.
#' @param negData A logical value indicating whether the data has both
#' positive and negative scores and the y-axis should be split.
#' Default value is \code{negData = FALSE}.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param linecolor A character value or vector of length 2 specifying the
#' line color(s) outlining the signal track(s).
#' Default value is \code{linecolor = "#37a7db"}.
#' @param fill A character value or vector of length 2 specifying
#' the fill color(s) of the signal track(s). Default value is \code{fill = NA}.
#' @param ymax A numeric specifying the fraction of the max y-value
#' to set as the height of the plot. Default value is \code{ymax = 1}.
#' @param range A numeric vector of length 2 specifying the y-range
#' of data to plot (c(min, max)).
#' @param scale A logical value indicating whether to include a data
#' scale label in the top left corner of the plot.
#' Default value is \code{scale = FALSE}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a
#' baseline along the x-axis. Default value is \code{baseline = TRUE}.
#' @param baseline.color Baseline color.
#' Default value is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width.
#' Default value is \code{baseline.lwd = 1}.
#' @param orientation A string specifying signal track orientation.
#' Default value is \code{orientation = "h"}. Options are:
#' \itemize{
#' \item{\code{"v"}: }{Vertical signal track orientation.}
#' \item{\code{"h"}: }{Horizontal signal track orientation.}
#' }
#' @param x A numeric or unit object specifying signal plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying signal plot y-location.
#' The character value will
#' place the signal plot y relative to the bottom of the most recently
#' plotted plot according to the units of the plotgardener page.
#' @param width A numeric or unit object specifying signal plot width.
#' @param height A numeric or unit object specifying signal plot height.
#' @param just Justification of signal plot relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal justification
#' and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if
#' \code{x}, \code{y}, \code{width}, or \code{height} are only given as
#' numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be
#' produced. Default value \code{draw = TRUE}.
#' @param params An optional \link[plotgardener]{pgParams} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{signal} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load signal data
#' library(plotgardenerData)
#' data("IMR90_ChIP_H3K27ac_signal")
#' data("GM12878_ChIP_H3K27ac_signal")
#'
#' ## Create a page
#' pageCreate(width = 7.5, height = 2.1, default.units = "inches")
#'
#' ## Define region
#' region <- pgParams(
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     range = c(0, 45)
#' )
#'
#' ## Plot and place signal plots
#' signal1 <- plotSignal(
#'     data = IMR90_ChIP_H3K27ac_signal, params = region,
#'     x = 0.5, y = 0.25, width = 6.5, height = 0.65,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' signal2 <- plotSignal(
#'     data = GM12878_ChIP_H3K27ac_signal, params = region,
#'     linecolor = "#7ecdbb",
#'     x = 0.5, y = 1, width = 6.5, height = 0.65,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Plot genome label
#' plotGenomeLabel(
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     x = 0.5, y = 1.68, length = 6.5,
#'     default.units = "inches"
#' )
#'
#' ## Add text labels
#' plotText(
#'     label = "IMR90", fonsize = 10, fontcolor = "#37a7db",
#'     x = 0.5, y = 0.25, just = c("left", "top"),
#'     default.units = "inches"
#' )
#' plotText(
#'     label = "GM12878", fonsize = 10, fontcolor = "#7ecdbb",
#'     x = 0.5, y = 1, just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' #A signal track can be placed on a plotgardener coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' plotSignal(data, chrom,
#'             chromstart = NULL, chromend = NULL,
#'             x, y, width, height, just = c("left", "top"),
#'             default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated
#' signal track by ignoring plot placement parameters:
#' \preformatted{
#' plotSignal(data, chrom,
#'             chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
plotSignal <- function(data, binSize = NA, binCap = TRUE, negData = FALSE,
                        chrom, chromstart = NULL, chromend = NULL,
                        assembly = "hg38", linecolor = "#37a7db",
                        fill = NA, ymax = 1, range = NULL, scale = FALSE,
                        bg = NA, baseline = TRUE, baseline.color = "grey",
                        baseline.lwd = 1, orientation = "h",
                        x = NULL, y = NULL, width = NULL,
                        height = NULL, just = c("left", "top"),
                        default.units = "inches", draw = TRUE,
                        params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for plotSignal
    errorcheck_plotSignal <- function(signal, signaltrack, fill) {
        dfChecks <- function(signal) {
            if (!"data.frame" %in% class(signal)) {
                if (!"GRanges" %in% class(signal)) {
                    if (!file_ext(signal) %in% c(
                        "bw", "bigWig",
                        "bigwig", "bedgraph"
                    )) {
                        stop("Invalid input. File must have a valid bigwig or ",
                            "bedgraph extension.", call. = FALSE)
                    }
                }
            }
        }

        if (is(signal, "list")) {
            if (length(signal) > 2) {
                stop("Invalid signal input. More than 2 signals provided.",
                    call. = FALSE
                )
            }

            invisible(lapply(signal, dfChecks))
        } else {
            dfChecks(signal = signal)
        }

        ## Genomic region
        regionErrors(chromstart = signaltrack$chromstart,
                        chromend = signaltrack$chromend)

        ## Range errors
        rangeErrors(range = signaltrack$range)
        
        checkColorby(fill = fill,
                        colorby = FALSE)
    }

    ## Define a function that reads in signal data for plotSignal
    read_signal <- function(signal, signaltrack) {

        signal <- read_rangeData(data = signal,
                                assembly = signaltrack$assembly,
                                chrom = signaltrack$chrom,
                                start = signaltrack$chromstart,
                                end = signaltrack$chromend)
        
        ## Check for overlapping data ranges
        if (any(IRanges::overlapsAny(
            GenomicRanges::makeGRangesFromDataFrame(signal),
            drop.self = TRUE
        ) == TRUE)) {
            stop("Data ranges cannot overlap. Please check `start` ",
                "and `end` column ranges.", call. = FALSE)
        }

        return(signal)
    }

    ## Define a function that formats/filters signal data
    format_data <- function(signal, signaltrack) {
        if (!any(colnames(signal) == "score")) {
            stop("Cannot find associated `score` column in data.",
                call. = FALSE
            )
        } else {
            ## Grab chrom, start, end and a "score" column
            signal <- signal[c("chrom", "start", "end", "score")]
        }

        ## Ensure the chromosome is a character
        signal[, "chrom"] <- as.character(signal[, "chrom"])

        ## Filter for desired region
        signal <- signal[
            which(signal[, "chrom"] == signaltrack$chrom &
                ((signal[, "start"] > signaltrack$chromstart &
                    signal[, "start"] < signaltrack$chromend |
                    signal[, "end"] > signaltrack$chromstart &
                        signal[, "end"] < signaltrack$chromend))),
            (2:4)
        ]
        ## Remove any duplicate rows
        signal <- signal[!duplicated(signal), ]

        ## Remove any NaN score values
        signal <- na.omit(signal)

        return(signal)
    }

    ## Define a function that checks and adjust the number/sizes of bins
    check_binNum <- function(signaltrack, binCap) {
        if (!is.na(signaltrack$binSize)) {
            if (signaltrack$binSize %% 0.25 != 0) {
                updated_binSize <- round(signaltrack$binSize / 0.25) * 0.25
                signaltrack$binSize <- updated_binSize
            }


            binNum <- (signaltrack$chromend - signaltrack$chromstart) /
                signaltrack$binSize
            signaltrack$binNum <- binNum

            if (!is.nan(binNum)) {

                ## Scale back binNum and print warning if binNum is
                ## greater than 8000
                if (binNum > 8000 && binCap == TRUE) {
                    updated_binNum <- 8000
                    updated_binSize <- (signaltrack$chromend -
                        signaltrack$chromstart) / binNum
                    signaltrack$binNum <- updated_binNum
                    signaltrack$binSize <- updated_binSize
                    warning("Too many bins: adjusting to 8000 bins of size ",
                        binSize, ". To override try binCap = FALSE.",
                        call. = FALSE
                    )
                }

                ## Scale bin size to 1 if binNum is larger than span
                if (binNum > (signaltrack$chromend - signaltrack$chromstart)) {
                    updated_binNum <- (signaltrack$chromend -
                        signaltrack$chromstart)
                    updated_binSize <- 1
                    signaltrack$binNum <- updated_binNum
                    signaltrack$binSize <- updated_binSize
                    warning("Number of bins larger than plot length: ",
                        "adjusting to ", binNum, " bins of size 1.",
                        call. = FALSE
                    )
                }
            } else {
                signaltrack$binSize <- NA
            }
        }

        return(signaltrack)
    }

    ## Define a function that bins, links, sorts, and combines data
    parseData <- function(signal, signaltrack) {
        if (!is.na(signaltrack$binSize)) {

            # =================================================================
            # BIN DATA
            # =================================================================
            ## Find the max signal value for each bin

            binChromend <- signaltrack$binSize * signaltrack$binNum +
                signaltrack$chromstart + signaltrack$binSize

            binDF <- data.frame(
                "start" = seq(
                    signaltrack$chromstart,
                    binChromend - signaltrack$binSize,
                    signaltrack$binSize
                ),
                "end" = seq(
                    signaltrack$chromstart + signaltrack$binSize,
                    binChromend,
                    signaltrack$binSize
                )
            )
            binDF$score <- rebinBigwig(signal, binDF)

            ## Use binned data as new signal data
            newSignal <- binDF
            # =================================================================
            # LINKING REGIONS
            # =================================================================

            linking_regions <- cbind(
                "start" = newSignal[seq(1, (nrow(newSignal) - 1)), "end"],
                "end" = newSignal[seq(2, nrow(newSignal)), "start"]
            )
            
            linking_regions <- matrix(linking_regions[which(
                linking_regions[, "start"] != linking_regions[, "end"]
            ), ], ncol = 2, dimnames = list(NULL, c("start", "end")))

            if (nrow(linking_regions) > 0) {
                linking_regions <- cbind(linking_regions, "score" = 0)
                
                ## Add linking regions to signaltrack
                newSignal <- rbind(newSignal, linking_regions)
            }

            # =================================================================
            # SORT AND COMBINE DATA
            # =================================================================

            ## Sort data
            newSignal <- newSignal[order(newSignal[, "start"]), ]
            
            ## Convert two columns to one
            newSignal <- cbind(
                "x" = as.vector(t(newSignal[, c("start", "end")])),
                "score" = as.vector(t(newSignal[, c("score", "score")]))
            )

            signal <- newSignal
        }


        return(signal)
    }

    ## Define a function that adjusts the range
    set_range <- function(signal1, signal2, signaltrack, split, pos = TRUE) {
        if (split == TRUE) {

            ## posSignal
            if (pos == TRUE) {
                if (is.null(signaltrack$range)) {
                    if (nrow(signal1) >= 2 & max(signal2[,"score"]) > 0) {
                        signaltrack$range[2] <- signaltrack$ymax *
                            max(signal2[, "score"])
                    } else {
                        signaltrack$range[2] <- 1
                    }
                }
            }

            ## negSignal
            if (pos == FALSE) {
                if (is.na(signaltrack$range[1])) {
                    if (nrow(signal1) >= 2 & min(signal2[,"score"]) < 0) {
                        signaltrack$range[1] <- signaltrack$ymax *
                            min(signal2[, "score"])
                    } else {
                        signaltrack$range[1] <- -1
                    }
                }
            }
        } else {
            if (is.null(signaltrack$range)) {
                if (nrow(signal1) >= 2 & max(signal2[,"score"]) > 0) {
                    
                    signaltrack$range <- c(0, signaltrack$ymax *
                        max(signal2[, "score"]))
                } else {
                    signaltrack$range <- c(0, 1)
                }
            }
        }

        return(signaltrack)
    }

    ## Define a function that parses out one vs. two fillcolors/linecolors
    parseColors <- function(color) {
        if (length(color) >= 2) {
            posCol <- color[1]
            negCol <- color[2]
        } else {
            posCol <- color
            negCol <- color
        }

        return(list(posCol, negCol))
    }

    ## Define a function that makes a grob for a signal (pos/neg)
    sigGrob <- function(signal, fillCol, lineCol, gp) {
        gp$col <- lineCol
        if (!is.null(fillCol) & !is.na(fillCol)) {
            if ("alpha" %in% names(gp)) {
                fillCol <- makeTransparent(color = fillCol, alpha = gp$alpha)
            }

            gp$fill <- fillCol

            sigGrob <- polygonGrob(
                x = c(
                    signal[1, "x"], signal[, "x"],
                    signal[nrow(signal), "x"]
                ),
                y = c(0, signal[, "score"], 0), gp = gp,
                default.units = "native"
            )
        } else {
            sigGrob <- segmentsGrob(
                x0 = signal[seq(1, length(signal[, "x"]) - 1), "x"],
                y0 = signal[seq(1, length(signal[, "score"]) - 1), "score"],
                x1 = signal[seq(2, length(signal[, "x"])), "x"],
                y1 = signal[seq(2, length(signal[, "score"])), "score"],
                gp = gp, default.units = "native"
            )
        }


        ## Add grob to gtree
        assign("signal_grobs",
            addGrob(
                gTree = get("signal_grobs", envir = pgEnv),
                child = sigGrob
            ),
            envir = pgEnv
        )
    }

    ## Define a function that finds data that falls out of the plot's
    ## range and draws a tiny black line to indicate it
    cutoffGrobs <- function(signal, signaltrack, side) {
        grobCutoffs <- function(df, side) {
            x0 <- df[1]
            x1 <- df[2]

            if (side == "top") {
                y <- 1
            } else {
                y <- 0
            }

            cutoffGrob <- segmentsGrob(
                x0 = x0, x1 = x1, y0 = unit(y, "npc"),
                y1 = unit(y, "npc"),
                gp = gpar(lwd = 1, col = "grey"),
                default.units = "native"
            )
            assign("signal_grobs",
                addGrob(
                    gTree = get("signal_grobs", envir = pgEnv),
                    child = cutoffGrob
                ),
                envir = pgEnv
            )
        }

        if (side == "top") {
            outsideData <- which(signal[, "score"] > signaltrack$range[2])
        } else {
            outsideData <- which(signal[, "score"] < signaltrack$range[1])
        }

        ## Get index pairs of outside data
        outsidePairs <- as.integer(outsideData + 1)
        ## Combine and order
        Outsidei <- c(outsideData, outsidePairs)
        Outsidei <- Outsidei[order(Outsidei)]
        ## Get xcoords
        signal <- as.data.frame(signal)
        Outside <- signal[Outsidei, 1]
        ## x0s are odd indeces and x1s are even
        x0s <- Outside[c(TRUE, FALSE)]
        x1s <- Outside[c(FALSE, TRUE)]
        pairs <- data.frame("x0" = x0s, "x1" = x1s)

        invisible(apply(pairs, 1, grobCutoffs, side = side))
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    sigInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "sigInternal"
    )

    ## Set gp
    sigInternal$gp <- setGP(
        gpList = gpar(),
        params = sigInternal, ...
    )
    
    ## Justification
    sigInternal$just <- justConversion(just = sigInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    signal_track <- structure(list(
        chrom = sigInternal$chrom,
        chromstart = sigInternal$chromstart,
        chromend = sigInternal$chromend,
        assembly = sigInternal$assembly,
        binSize = sigInternal$binSize,
        binNum = NULL, range = sigInternal$range,
        ymax = sigInternal$ymax,
        x = sigInternal$x, y = sigInternal$y,
        width = sigInternal$width,
        height = sigInternal$height,
        just = sigInternal$just, grobs = NULL
    ),
    class = "signal"
    )
    attr(x = signal_track, which = "plotted") <- sigInternal$draw

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(sigInternal$data)) stop("argument \"data\" is missing, ",
                                        "with no default.", call. = FALSE)
    if (is.null(sigInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)

    check_placement(object = signal_track)
    errorcheck_plotSignal(
        signal = sigInternal$data,
        signaltrack = signal_track,
        fill = sigInternal$fill
    )

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    signal_track$assembly <- parseAssembly(assembly = signal_track$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    signal_track <- defaultUnits(
        object = signal_track,
        default.units = sigInternal$default.units
    )

    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================

    scaleChecks <- genomicScale(object = signal_track,
                                objectInternal = sigInternal,
                                plotType = "signal track")
    signal_track <- scaleChecks[[1]]
    sigInternal <- scaleChecks[[2]]

    # =========================================================================
    # SET BINSIZE
    # =========================================================================

    if (is.na(signal_track$binSize) == TRUE) {
        if (!is.null(signal_track$chromstart)) {
            binSize <- (signal_track$chromend - signal_track$chromstart) / 2000
            signal_track$binSize <- binSize
        }
    }

    # =========================================================================
    # CHECK AND ADJUST BIN NUMBER/BIN SIZE
    # =========================================================================

    signal_track <- check_binNum(
        signaltrack = signal_track,
        binCap = sigInternal$binCap
    )

    # =========================================================================
    # READ IN, FORMAT, FILTER, BIN, LINK AND SORT DATA
    # =========================================================================
    
    if (is(sigInternal$data, "list")) {
        signal <- lapply(sigInternal$data, read_signal,
            signaltrack = signal_track
        )
        signal <- lapply(signal, format_data, signaltrack = signal_track)
        posSignal <- signal[[1]]
        negSignal <- signal[[2]]
        split <- TRUE
    } else {
        signal <- read_signal(
            signal = sigInternal$data,
            signaltrack = signal_track
        )
        signal <- format_data(signal = signal, signaltrack = signal_track)

        if (any(signal[, "score"] < 0)) {
            if (sigInternal$negData == FALSE) {
                warning("Negative scores detected in signal data. To make ",
                "an entirely positive signal track, ",
                "please remove negative scores from data.", call. = FALSE)
            }

            posSignal <- signal[which(signal[, "score"] >= 0), ]
            negSignal <- signal[which(signal[, "score"] < 0), ]
            negSignal[, "score"] <- negSignal[, "score"] * -1
            split <- TRUE
        } else {
            posSignal <- signal
            split <- FALSE
            if (sigInternal$negData == TRUE) {
                negSignal <- data.frame()
                split <- TRUE
            }
        }
    }


    # =========================================================================
    # BIN, LINK, AND SORT DATA AND FIX Y-LIMITS
    # =========================================================================

    if (split == TRUE) {
        if (nrow(posSignal) >= 2) {
            posSignal2 <- parseData(
                signal = posSignal,
                signaltrack = signal_track
            )
        }

        signal_track <- set_range(
            signal1 = posSignal, signal2 = posSignal2,
            signaltrack = signal_track, split = TRUE
        )

        if (nrow(negSignal) >= 2) {
            negSignal2 <- parseData(
                signal = negSignal,
                signaltrack = signal_track
            )
            negSignal2[, "score"] <- negSignal2[, "score"] * -1
        }


        signal_track <- set_range(
            signal1 = negSignal, signal2 = negSignal2,
            signaltrack = signal_track, split = TRUE,
            pos = FALSE
        )


        if (signal_track$range[1] == 0 & signal_track$range[2] == 0) {
            signal_track$range <- c(-1, 1)
        }
    } else {
        if (nrow(posSignal) >= 2) {
            posSignal2 <- parseData(
                signal = posSignal,
                signaltrack = signal_track
            )
        }

        signal_track <- set_range(
            signal1 = posSignal, signal2 = posSignal2,
            signaltrack = signal_track, split = FALSE
        )
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "signal",
        length(grep(
            pattern = "signal",
            x = currentViewports
        )) + 1
    )

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(signal_track$x) | is.null(signal_track[["y"]])) {

        if (sigInternal$orientation == "h"){

            vp <- viewport(
                height = unit(0.25, "snpc"), width = unit(1, "snpc"),
                x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                clip = "on",
                xscale = sigInternal$xscale,
                yscale = c(signal_track$range[1], signal_track$range[2]),
                just = "center",
                name = paste0(vp_name, "_h")
                )
            } else if (sigInternal$orientation == "v"){

                ## outside clipping viewport
                vpClip <- viewport(
                    x = unit(0.5, "npc"),
                    y = unit(0.5, "npc"),
                    width = unit(0.25, "snpc"),
                    height = unit(1, "snpc"),
                    just = "center",
                    clip = "on",
                    xscale = c(signal_track$range[2], signal_track$range[1]),
                    yscale = sigInternal$xscale,
                    name = paste0(vp_name, "_vClip")

                )
                pushViewport(vpClip)
                height <- convertWidth(unit(1, "npc"), unitTo = "inches")
                width <- convertHeight(unit(1, "npc"), unitTo = "inches")
                upViewport()
                ## Make rotated, horizontal viewport
                vp <- viewport(
                    height = height, width = width,
                    x = unit(1, "npc"), y = unit(0, "npc"),
                    just = c("left", "bottom"),
                    xscale = sigInternal$xscale,
                    yscale = c(signal_track$range[1], signal_track$range[2]),
                    name = paste0(vp_name, "_v"),
                    angle = 90
                )
            }

        if (sigInternal$draw == TRUE) {
            vp$name <- "signal1"
            if (sigInternal$orientation == "v"){
                vpClip$name <- "signal1_Clip"
            }

            grid.newpage()
        }
        
        } else {

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = signal_track)

        if (sigInternal$orientation == "h"){

            ## Make viewport
            vp <- viewport(
                height = page_coords$height, width = page_coords$width,
                x = page_coords$x, y = page_coords$y,
                clip = "on",
                xscale = sigInternal$xscale,
                yscale = c(signal_track$range[1], signal_track$range[2]),
                just = sigInternal$just,
                name = paste0(vp_name, "_h")
            )
            addViewport(paste0(vp_name, "_h"))
        } else if (sigInternal$orientation == "v"){

            ## outside clipping viewport
            vpClip <- viewport(
                x = page_coords$x,
                y = page_coords$y,
                width = page_coords$width,
                height = page_coords$height,
                just = sigInternal$just,
                clip = "on",
                xscale = c(signal_track$range[2], signal_track$range[1]),
                yscale = sigInternal$xscale,
                name = paste0(vp_name, "_vClip")

            )
            ## Make rotated, horizontal viewport
            vp <- viewport(
                height = page_coords$width, width = page_coords$height,
                x = unit(1, "npc"), y = unit(0, "npc"),
                just = c("left", "bottom"),
                xscale = sigInternal$xscale,
                yscale = c(signal_track$range[1], signal_track$range[2]),
                name = paste0(vp_name, "_v"),
                angle = 90
            )
            addViewport(paste0(vp_name, "_vClip"))
        }

    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================

    backgroundGrob <- rectGrob(gp = gpar(
        fill = sigInternal$bg,
        col = NA
    ), name = "background")
    assign("signal_grobs", gTree(vp = vp, children = gList(backgroundGrob)),
        envir = pgEnv
    )

    # =========================================================================
    # MAKE GROBS
    # =========================================================================

    if (!is.na(signal_track$binSize)) {
        if (split == TRUE) {
            fills <- parseColors(sigInternal$fill)
            lines <- parseColors(sigInternal$linecolor)

            if (nrow(posSignal) >= 2) {
                sigGrob(
                    signal = posSignal2, fillCol = fills[[1]],
                    lineCol = lines[[1]], gp = sigInternal$gp
                )
                ## Find and make cutoff lines
                cutoffGrobs(
                    signal = posSignal2, signaltrack = signal_track,
                    side = "top"
                )
            } else {
                sigInternal$gp$col <- lines[[1]]
                posGrob <- segmentsGrob(
                    x0 = 0, y0 = unit(0, "native"),
                    x1 = 1, y1 = unit(0, "native"),
                    gp = sigInternal$gp
                )
                assign("signal_grobs",
                    addGrob(
                        gTree = get("signal_grobs", envir = pgEnv),
                        child = posGrob
                    ),
                    envir = pgEnv
                )
                warning("Not enough top signal data to plot.", call. = FALSE)
            }


            if (nrow(negSignal) >= 2) {
                sigGrob(
                    signal = negSignal2, fillCol = fills[[2]],
                    lineCol = lines[[2]], gp = sigInternal$gp
                )
                ## Find and make cutoff lines
                cutoffGrobs(
                    signal = negSignal2, signaltrack = signal_track,
                    side = "bottom"
                )
            } else {
                sigInternal$gp$col <- lines[[2]]
                negGrob <- segmentsGrob(
                    x0 = 0, y0 = unit(0, "native"),
                    x1 = 1, y1 = unit(0, "native"),
                    gp = sigInternal$gp
                )
                assign("signal_grobs",
                    addGrob(
                        gTree = get("signal_grobs", envir = pgEnv),
                        child = negGrob
                    ),
                    envir = pgEnv
                )
                warning("Not enough bottom signal data to plot.",
                    call. = FALSE
                )
            }


            lineGrob <- segmentsGrob(
                x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                y0 = 0, y1 = 0,
                gp = gpar(
                    col = sigInternal$baseline.color,
                    lwd = sigInternal$baseline.lwd
                ),
                default.units = "native"
            )
            assign("signal_grobs",
                addGrob(
                    gTree = get("signal_grobs", envir = pgEnv),
                    child = lineGrob
                ),
                envir = pgEnv
            )
        } else {
            if (nrow(posSignal) >= 2) {
                if (sigInternal$baseline == TRUE) {
                    baselineGrob <- segmentsGrob(
                        x0 = unit(0, "npc"),
                        x1 = unit(1, "npc"),
                        y0 = 0, y1 = 0,
                        gp = gpar(
                            col = sigInternal$baseline.color,
                            lwd = sigInternal$baseline.lwd
                        ),
                        default.units = "native"
                    )
                    assign("signal_grobs",
                        addGrob(
                            gTree = get("signal_grobs", envir = pgEnv),
                            child = baselineGrob
                        ),
                        envir = pgEnv
                    )
                }

                sigGrob(
                    signal = posSignal2, fillCol = sigInternal$fill[1],
                    lineCol = sigInternal$linecolor[1],
                    gp = sigInternal$gp
                )
                ## Find and make cutoff lines
                cutoffGrobs(
                    signal = posSignal2, signaltrack = signal_track,
                    side = "top"
                )
            } else {
                sigInternal$gp$col <- sigInternal$linecolor
                signalGrob <- segmentsGrob(
                    x0 = 0, y0 = unit(0, "native"),
                    x1 = 1, y1 = unit(0, "native"),
                    gp = sigInternal$gp
                )
                assign("signal_grobs",
                    addGrob(
                        gTree = get("signal_grobs", envir = pgEnv),
                        child = signalGrob
                    ),
                    envir = pgEnv
                )
                warning("Not enough data within range to plot.", call. = FALSE)
            }
        }
        # =====================================================================
        # SCALE
        # =====================================================================

        ## Add scale of the range of data in the top left corner
        if (sigInternal$scale == TRUE) {
            upperLim <- round(signal_track$range[2], digits = 4)
            lowerLim <- round(signal_track$range[1], digits = 4)
            scaleGrob <- textGrob(
                label = paste0("[", lowerLim, " - ", upperLim, "]"),
                just = c("left", "top"), x = 0, y = 1,
                gp = sigInternal$gp
            )

            ## Add grob to gtree
            assign("signal_grobs",
                addGrob(
                    gTree = get("signal_grobs", envir = pgEnv),
                    child = scaleGrob
                ),
                envir = pgEnv
            )
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (sigInternal$draw == TRUE) {

        if (sigInternal$orientation == "v"){
            pushViewport(vpClip)
            grid.draw(get("signal_grobs", envir = pgEnv))
            upViewport()
        } else if (sigInternal$orientation == "h"){
            grid.draw(get("signal_grobs", envir = pgEnv))
        }

    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    signal_track$grobs <- get("signal_grobs", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("signal[", vp_name, "]")
    invisible(signal_track)
}
