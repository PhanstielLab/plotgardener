#' Plot paired-end genomic range data in an arch style
#' 
#' @usage bbPlotPairsArches(
#'     data,
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
#'     style = "2D",
#'     flip = FALSE,
#'     curvature = 5,
#'     archHeight = NULL,
#'     fill = "#1f4297",
#'     linecolor = NA,
#'     alpha = 0.4,
#'     bg = NA,
#'     clip = FALSE,
#'     baseline = FALSE,
#'     baseline.color = "grey",
#'     baseline.lwd = 1,
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
#' @param data A string specifying the BEDPE file path, a dataframe
#' in BEDPE format specifying data to be plotted, or a
#' \link[InteractionSet]{GInteractions} object.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[BentoBox]{bbAssembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param style Character value describing the style of arches.
#' Default value is \code{style = "2D"}. Options are:
#' \itemize{
#' \item{\code{"2D"}: }{Arches will be drawn in a 2-dimensional style.}
#' \item{\code{"3D"}: }{Arches will be drawn in a 3-dimensional style.}
#' }
#' @param flip Logical value indicating whether to reflect arches over
#' the x-axis. Default value is \code{flip = FALSE}.
#' @param curvature Numeric indicating the number of points along the
#' arch curvature. Default value is \code{curvature = 5}.
#' @param archHeight Single numeric value or numeric vector specifying
#' the arch heights. When NULL, all arches will be the same height,
#' filling up the given plot area
#' @param fill A single character value, a vector, or a 
#' \link[BentoBox]{colorby} object specifying fill colors of arches.
#' Default value is \code{fill = #1f4297"}.
#' @param linecolor A single character value, a vector, or a
#' \link[BentoBox]{colorby} object specifying the color of the lines
#' outlining arches. Default value is \code{linecolor = NA}.
#' Special options include:
#' \itemize{
#' \item{\code{NA}: }{No line color.}
#' \item{\code{"fill"}: }{Same color as \code{fill}.}
#' }
#' @param alpha Numeric value specifying transparency.
#' Default value is \code{alpha = 0.4}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param clip A logical value indicating whether to clip any
#' arches that get cutoff in the given genomic region.
#' Default value is \code{clip = FALSE}.
#' @param baseline Logical value indicating whether to include
#' a baseline along the x-axis. Default value is \code{baseline = FALSE}.
#' @param baseline.color Baseline color.
#' Default value is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width.
#' Default value is \code{baseline.lwd = 1}.
#' @param x A numeric or unit object specifying pair arches plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying BEDPE arches plot y-location.
#' The character value will
#' place the pair arches plot y relative to the bottom of the most
#' recently plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying pair arches plot width.
#' @param height A numeric or unit object specifying pair arches plot height.
#' @param just Justification of pair arches plot relative to its (x, y)
#' location. If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if
#' \code{x}, \code{y}, \code{width}, or \code{height} are only given as
#' numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be
#' produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bbParams} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_arches} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load paired ranges data in BEDPE format
#' library(BentoBoxData)
#' data("IMR90_DNAloops_pairs")
#'
#' ## Set the coordinates
#' params <- bbParams(
#'     chrom = "chr21",
#'     chromstart = 27900000, chromend = 30700000,
#'     assembly = "hg19",
#'     width = 7
#' )
#'
#' ## Create a page
#' bbPageCreate(width = 7.5, height = 2.1, default.units = "inches")
#'
#' ## Add a length column to color by
#' IMR90_DNAloops_pairs$length <- 
#'         (IMR90_DNAloops_pairs$start2 - IMR90_DNAloops_pairs$start1) / 1000
#'
#' ## Translate lengths into heights
#' heights <- IMR90_DNAloops_pairs$length / max(IMR90_DNAloops_pairs$length)
#'
#' ## Plot the data
#' archPlot <- bbPlotPairsArches(
#'     data = IMR90_DNAloops_pairs, params = params,
#'     fill = colorby("length", palette = 
#'                 colorRampPalette(c("dodgerblue2", "firebrick2"))),
#'     linecolor = "fill",
#'     archHeight = heights, alpha = 1,
#'     x = 0.25, y = 0.25, height = 1.5,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bbAnnoGenomeLabel(plot = archPlot, x = 0.25, y = 1.78, scale = "Mb")
#'
#'
#' ## Annotate heatmap legend
#' bbAnnoHeatmapLegend(
#'     plot = archPlot, fontcolor = "black",
#'     x = 7.0, y = 0.25,
#'     width = 0.10, height = 1, fontsize = 10
#' )
#'
#' ## Add the heatmap legend title
#' bbPlotText(
#'     label = "Kb", rot = 90, x = 6.9, y = 0.75,
#'     just = c("center", "center"),
#'     fontsize = 10
#' )
#'
#'
#' ## Hide page guides
#' bbPageGuideHide()
#' @details
#' A pair arches plot can be placed on a BentoBox coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' bbPlotPairsArches(data chrom,
#'                 chromstart = NULL, chromend = NULL,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated pair
#' arches plot by ignoring plot placement parameters:
#' \preformatted{
#' bbPlotPairsArches(data, chrom,
#'                 chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
bbPlotPairsArches <- function(data, chrom, chromstart = NULL, chromend = NULL,
                            assembly = "hg38", style = "2D", flip = FALSE,
                            curvature = 5, archHeight = NULL,
                            fill = "#1f4297",
                            linecolor = NA, alpha = 0.4, bg = NA,
                            clip = FALSE, baseline = FALSE,
                            baseline.color = "grey", baseline.lwd = 1,
                            x = NULL, y = NULL, width = NULL, height = NULL,
                            just = c("left", "top"),
                            default.units = "inches", draw = TRUE,
                            params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors
    errorcheck_bbArches <- function(bedpe, archesPlot, style, fill) {
        ## Genomic region
        bb_regionErrors(chromstart = archesPlot$chromstart,
                    chromend = archesPlot$chromend)

        if (!style %in% c("3D", "2D")) {
            stop("Invalid \'style\' input. Options are \'3D\' and \'2D\'.",
                call. = FALSE
            )
        }

        ## Fill colorby checks
        bb_checkColorby(fill = fill,
                        colorby = TRUE,
                        data = bedpe)
    }

    ## Define a function that will produce a yscale for arches based on height
    height_yscale <- function(heights, flip) {
        if (length(heights) == 1) {
            if (flip == FALSE) {
                yscale <- c(0, heights)
            } else {
                yscale <- c(heights, 0)
            }
        } else {
            if (flip == FALSE) {
                yscale <- c(0, max(heights))
            } else {
                yscale <- c(max(heights), 0)
            }
        }
        return(yscale)
    }

    ## Define a function that normalizes arch heights
    normHeights <- function(height, min, max) {

        ## First normalize from 0 to 1
        newHeight <- (height - min) / (max - min)

        ## Then scale to a range of 1.38
        finalHeight <- newHeight * 1.38

        return(finalHeight)
    }

    ## Define a function that creates ribbon arch grobs
    drawRibbons <- function(df, style, arch, flip, transp, gp) {
        x1 <- df$start1
        x2 <- df$end1
        y1 <- df$start2
        y2 <- df$end2
        fillCol <- df$color
        lineCol <- df$linecolor
        outerHeight <- df$normHeight
        innerHeight <- outerHeight - 0.01
        gp$fill <- fillCol
        gp$col <- lineCol
        gp$alpha <- transp

        if (style == "3D") {
            x1 <- df$end1
            x2 <- df$start1
        }

        ## Designate bezier control points

        innerX <- unit(
            seq(x2, y1, length.out = arch)[c(1, 2, seq((arch - 1), arch))],
            "native"
        )
        outerX <- unit(
            seq(x1, y2, length.out = arch)[c(1, 2, seq((arch - 1), arch))],
            "native"
        )

        if (flip == FALSE) {
            ## Switch y-positions for top plotting
            innerY <- unit(c(0, innerHeight, innerHeight, 0), "npc")
            outerY <- unit(c(0, outerHeight, outerHeight, 0), "npc")
        } else {
            ## Switch y-positions for bottom plotting
            innerY <- unit(c(1, 1 - innerHeight, 1 - innerHeight, 1), "npc")
            outerY <- unit(c(1, 1 - outerHeight, 1 - outerHeight, 1), "npc")
        }

        ## Calculate loop arcs using bezier curves
        innerLoop <- bezierGrob(x = innerX, y = innerY)
        outerLoop <- bezierGrob(x = outerX, y = outerY)

        ## Extract points from bezier curves
        innerBP <- bezierPoints(innerLoop)
        outerBP <- bezierPoints(outerLoop)

        ## Connect points, convert to proper units and draw polygons
        archGrob <- polygonGrob(
            x = unit(
                c(
                    convertX(outerBP$x, "native"),
                    rev(convertX(innerBP$x, "native"))
                ),
                "native"
            ),
            y = unit(
                c(
                    convertY(outerBP$y, "npc"),
                    rev(convertY(innerBP$y, "npc"))
                ),
                "npc"
            ),
            gp = gp
        )
        assign("arches_grobs",
            addGrob(
                gTree = get("arches_grobs", envir = bbEnv),
                child = archGrob
            ),
            envir = bbEnv
        )
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_archInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_archInternal"
    )

    ## Set gp
    bb_archInternal$gp <- setGP(
        gpList = gpar(),
        params = bb_archInternal, ...
    )
    
    ## Justification
    bb_archInternal$just <- bb_justConversion(just = bb_archInternal$just)

    # =========================================================================
    # CHECK ARGUMENT ERRORS
    # =========================================================================
    if (is.null(bb_archInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(bb_archInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    archesPlot <- structure(list(
        bedpe = NULL, chrom = bb_archInternal$chrom,
        chromstart = bb_archInternal$chromstart,
        chromend = bb_archInternal$chromend,
        assembly = bb_archInternal$assembly,
        color_palette = NULL,
        zrange = NULL,
        x = bb_archInternal$x, y = bb_archInternal$y,
        width = bb_archInternal$width,
        height = bb_archInternal$height,
        just = bb_archInternal$just, grobs = NULL
    ),
    class = "bb_arches"
    )
    attr(x = archesPlot, which = "plotted") <- bb_archInternal$draw

    # =========================================================================
    # CHECK PLACEMENT ERRORS
    # =========================================================================

    check_placement(object = archesPlot)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    archesPlot$assembly <- parse_bbAssembly(assembly = archesPlot$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    archesPlot <- defaultUnits(
        object = archesPlot,
        default.units = bb_archInternal$default.units
    )

    # =========================================================================
    # READ IN FILE OR DATAFRAME
    # =========================================================================

    bedpe <- read_pairedData(data = bb_archInternal$data,
                            assembly = archesPlot$assembly)
    
    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    errorcheck_bbArches(
        bedpe = bedpe, archesPlot = archesPlot,
        style = bb_archInternal$style,
        fill = bb_archInternal$fill
    )
    
    ## chrom format and data chrom format
    chromDataAgreement(data = bedpe, chrom = archesPlot$chrom,
                    type = "pairs")

    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================

    scaleChecks <- genomicScale(object = archesPlot,
                                objectInternal = bb_archInternal,
                                plotType = "paired data arches")
    archesPlot <- scaleChecks[[1]]
    bb_archInternal <- scaleChecks[[2]]
    
    # =========================================================================
    # COLORS
    # =========================================================================
    
    if (bb_archInternal$clip == TRUE){
        subset <- "pairs_clip"
    } else {
        subset <- "pairs"
    }
    
    archColors <- bb_parseColors(data = bedpe,
                                fill = bb_archInternal$fill,
                                object = archesPlot,
                                subset = subset)
    if (length(archColors[[1]]) > 0){
        bedpe$color <- archColors[[1]]
    } else {
        bedpe$color <- rep("#1f4297", nrow(bedpe))
    }
    
    archesPlot <- archColors[[2]]
    
    bedpe$linecolor <- bb_lineColors(linecolor = bb_archInternal$linecolor,
                                    fillcolors = bedpe$color,
                                    data = bedpe,
                                    object = archesPlot,
                                    subset = subset)
    
    # =========================================================================
    # SUBSET DATA
    # =========================================================================

    if (!is.null(archesPlot$chromstart) & !is.null(archesPlot$chromend)) {
        if (bb_archInternal$clip == TRUE) {
            bedpe <- bedpe[which(bedpe[, "chrom1"] == archesPlot$chrom &
                bedpe[, "chrom2"] == archesPlot$chrom &
                bedpe[, "start1"] >= archesPlot$chromstart &
                bedpe[, "end1"] <= archesPlot$chromend &
                bedpe[, "start2"] >= archesPlot$chromstart &
                bedpe[, "end2"] <= archesPlot$chromend), ]
        } else {
            bedpe <- bedpe[which(bedpe[, "chrom1"] == archesPlot$chrom &
                bedpe[, "chrom2"] == archesPlot$chrom &
                ((bedpe[, "end1"] >= archesPlot$chromstart &
                    bedpe[, "end1"] <= archesPlot$chromend) |
                    (bedpe[, "start2"] <= archesPlot$chromstart &
                        bedpe[, "start2"] >= archesPlot$chromend))), ]
        }
    } else {
        bedpe <- data.frame(matrix(nrow = 0, ncol = 6))
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_arches",
        length(grep(
            pattern = "bb_arches",
            x = currentViewports
        )) + 1
    )

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(archesPlot$x) | is.null(archesPlot$y)) {
        vp <- viewport(
            height = unit(0.5, "npc"), width = unit(1, "npc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            xscale = bb_archInternal$xscale,
            clip = "on",
            just = "center",
            name = vp_name
        )

        if (bb_archInternal$draw == TRUE) {
            vp$name <- "bb_arches1"
            grid.newpage()
        }
    } else {
        add_bbViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = archesPlot)

        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            xscale = bb_archInternal$xscale,
            clip = "on",
            just = bb_archInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # HEIGHTS
    # =========================================================================

    if (nrow(bedpe) > 0) {
        if (is.null(bb_archInternal$archHeight)) {
            bedpe$height <- rep(1, nrow(bedpe))
        } else if (length(bb_archInternal$archHeight) == 1) {
            bedpe$height <- rep(bb_archInternal$archHeight, nrow(bedpe))
            yscale <- height_yscale(
                heights = bb_archInternal$archHeight,
                flip = bb_archInternal$flip
            )
            vp$yscale <- yscale
        } else {
            bedpe$height <- bb_archInternal$archHeight[seq(1, nrow(bedpe))]
            yscale <- height_yscale(
                heights = bb_archInternal$archHeight,
                flip = bb_archInternal$flip
            )
            vp$yscale <- yscale
        }

        if (length(bedpe$height) > 0) {
            bedpe$normHeight <- lapply(bedpe$height, normHeights,
                min = 0,
                max = max(bedpe$height)
            )
        }
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================

    backgroundGrob <- rectGrob(
        gp = gpar(fill = bb_archInternal$bg, col = NA),
        name = "background"
    )
    assign("arches_grobs", gTree(vp = vp, children = gList(backgroundGrob)),
        envir = bbEnv
    )

    # =========================================================================
    # GROBS
    # =========================================================================

    if (nrow(bedpe) > 0) {
        if (bb_archInternal$baseline == TRUE) {
            baselineGrob <- segmentsGrob(
                x0 = unit(0, "npc"), y0 = 0,
                x1 = unit(1, "npc"), y1 = 0,
                default.units = "native",
                gp = gpar(
                    col = bb_archInternal$baseline.color,
                    lwd = bb_archInternal$baseline.lwd
                )
            )
            assign("arches_grobs",
                addGrob(
                    gTree = get("arches_grobs", envir = bbEnv),
                    child = baselineGrob
                ),
                envir = bbEnv
            )
        }
        
        invisible(apply(bedpe, 1, drawRibbons,
            style = bb_archInternal$style,
            arch = bb_archInternal$curvature,
            flip = bb_archInternal$flip,
            transp = bb_archInternal$alpha, gp = bb_archInternal$gp
        ))
    } else {
        if (bb_archInternal$txdbChecks == TRUE) {
            warning("Data contains no values.", call. = FALSE)
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (bb_archInternal$draw == TRUE) {
        grid.draw(get("arches_grobs", envir = bbEnv))
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    archesPlot$grobs <- get("arches_grobs", envir = bbEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_arches[", vp$name, "]")
    invisible(archesPlot)
}
