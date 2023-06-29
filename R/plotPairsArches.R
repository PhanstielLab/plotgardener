#' Plot paired-end genomic range data in an arch style
#' 
#' @usage plotPairsArches(
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
#'     clip.noAnchor = TRUE,
#'     range = NULL,
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
#' \link[plotgardener]{assembly} object.
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
#' @param archHeight Single numeric value, numeric vector, or column name 
#' in data specifying the arch heights. When NULL, all arches will be the 
#' same height, filling up the given plot area.
#' @param fill A single character value, a vector, or a 
#' \link[plotgardener]{colorby} object specifying fill colors of arches.
#' Default value is \code{fill = #1f4297"}.
#' @param linecolor A single character value, a vector, or a
#' \link[plotgardener]{colorby} object specifying the color of the lines
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
#' @param clip.noAnchor A logical value indicating whether to clip
#' any arches that overlap the given genomic region but do not 
#' have an anchor in that region. Default value is \code{clip.noAnchor = TRUE}.
#' @param range A numeric vector of length 2 specifying the y-range
#' of \code{archHeight} to plot (c(min, max)).
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
#' recently plotted plot according to the units of the plotgardener page.
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
#' @param params An optional \link[plotgardener]{pgParams} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{arches} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load paired ranges data in BEDPE format
#' library(plotgardenerData)
#' data("IMR90_DNAloops_pairs")
#'
#' ## Set the coordinates
#' params <- pgParams(
#'     chrom = "chr21",
#'     chromstart = 27900000, chromend = 30700000,
#'     assembly = "hg19",
#'     width = 7
#' )
#'
#' ## Create a page
#' pageCreate(width = 7.5, height = 2.1, default.units = "inches")
#'
#' ## Add a length column to color by
#' IMR90_DNAloops_pairs$length <- 
#'         (IMR90_DNAloops_pairs$start2 - IMR90_DNAloops_pairs$start1) / 1000
#'
#' ## Translate lengths into heights
#' IMR90_DNAloops_pairs$h <- 
#'         IMR90_DNAloops_pairs$length / max(IMR90_DNAloops_pairs$length)
#'
#' ## Plot the data
#' archPlot <- plotPairsArches(
#'     data = IMR90_DNAloops_pairs, params = params,
#'     fill = colorby("length", palette = 
#'                 colorRampPalette(c("dodgerblue2", "firebrick2"))),
#'     linecolor = "fill",
#'     archHeight = "h", alpha = 1,
#'     x = 0.25, y = 0.25, height = 1.5,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' annoGenomeLabel(plot = archPlot, x = 0.25, y = 1.78, scale = "Mb")
#'
#'
#' ## Annotate heatmap legend
#' annoHeatmapLegend(
#'     plot = archPlot, fontcolor = "black",
#'     x = 7.0, y = 0.25,
#'     width = 0.10, height = 1, fontsize = 10
#' )
#'
#' ## Add the heatmap legend title
#' plotText(
#'     label = "Kb", rot = 90, x = 6.9, y = 0.75,
#'     just = c("center", "center"),
#'     fontsize = 10
#' )
#'
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' A pair arches plot can be placed on a plotgardener coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' plotPairsArches(data chrom,
#'                 chromstart = NULL, chromend = NULL,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated pair
#' arches plot by ignoring plot placement parameters:
#' \preformatted{
#' plotPairsArches(data, chrom,
#'                 chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
plotPairsArches <- function(data, chrom, chromstart = NULL, chromend = NULL,
                            assembly = "hg38", style = "2D", flip = FALSE,
                            curvature = 5, archHeight = NULL,
                            fill = "#1f4297",
                            linecolor = NA, alpha = 0.4, bg = NA,
                            clip = FALSE, clip.noAnchor = TRUE, 
                            range = NULL, baseline = FALSE,
                            baseline.color = "grey", baseline.lwd = 1,
                            x = NULL, y = NULL, width = NULL, height = NULL,
                            just = c("left", "top"),
                            default.units = "inches", draw = TRUE,
                            params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors
    errorcheck_plotArches <- function(bedpe, archesPlot, style, fill) {
        ## Genomic region
        regionErrors(chromstart = archesPlot$chromstart,
                    chromend = archesPlot$chromend)

        if (!style %in% c("3D", "2D")) {
            stop("Invalid \'style\' input. Options are \'3D\' and \'2D\'.",
                call. = FALSE
            )
        }

        ## Fill colorby checks
        checkColorby(fill = fill,
                        colorby = TRUE,
                        data = bedpe)
        
        rangeErrors(range = archesPlot$range)
    }

    ## Define a function that will produce a yscale for arches based on height
    ## and range
    height_yscale <- function(heights, flip, range) {
        if (is.null(range)){
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
        } else {
            if (flip == FALSE){
                yscale <- range
            } else {
                yscale <- rev(range)
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

        x1 <- min(df$start1, df$start2)
        x2 <- min(df$end1, df$end2)
        y1 <- max(df$start1, df$start2)
        y2 <- max(df$end1, df$end2)
        
        
        fillCol <- df$color
        lineCol <- df$linecolor
        outerHeight <- df$normHeight
        gp$fill <- fillCol
        gp$col <- lineCol
        gp$alpha <- transp

        ## Calculate innerHeight
        anchor1 <- x2 - x1
        anchor2 <- y2 - y1
        if (anchor1 == anchor2){
            innerHeight <- outerHeight - 0.01
        } else {
            
            diff1 <- abs(y1-x2)
            diff2 <- abs(y2-x1)
            innerHeight <- outerHeight*((y1-x2)/(y2-x1))
        }

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
                gTree = get("arches_grobs", envir = pgEnv),
                child = archGrob
            ),
            envir = pgEnv
        )
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    archInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "archInternal"
    )

    ## Set gp
    archInternal$gp <- setGP(
        gpList = gpar(),
        params = archInternal, ...
    )
    
    ## Justification
    archInternal$just <- justConversion(just = archInternal$just)

    # =========================================================================
    # CHECK ARGUMENT ERRORS
    # =========================================================================
    if (is.null(archInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(archInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    archesPlot <- structure(list(
        bedpe = NULL, chrom = archInternal$chrom,
        chromstart = archInternal$chromstart,
        chromend = archInternal$chromend,
        assembly = archInternal$assembly,
        color_palette = NULL,
        zrange = NULL,
        range = archInternal$range,
        x = archInternal$x, y = archInternal$y,
        width = archInternal$width,
        height = archInternal$height,
        just = archInternal$just, grobs = NULL
    ),
    class = "arches"
    )
    attr(x = archesPlot, which = "plotted") <- archInternal$draw

    # =========================================================================
    # CHECK PLACEMENT ERRORS
    # =========================================================================

    check_placement(object = archesPlot)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    archesPlot$assembly <- parseAssembly(assembly = archesPlot$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    archesPlot <- defaultUnits(
        object = archesPlot,
        default.units = archInternal$default.units
    )

    # =========================================================================
    # READ IN FILE OR DATAFRAME
    # =========================================================================

    bedpe <- read_pairedData(data = archInternal$data,
                            assembly = archesPlot$assembly)
    
    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    errorcheck_plotArches(
        bedpe = bedpe, archesPlot = archesPlot,
        style = archInternal$style,
        fill = archInternal$fill
    )
    
    ## chrom format and data chrom format
    chromDataAgreement(data = bedpe, chrom = archesPlot$chrom,
                    type = "pairs")

    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================

    scaleChecks <- genomicScale(object = archesPlot,
                                objectInternal = archInternal,
                                plotType = "paired data arches")
    archesPlot <- scaleChecks[[1]]
    archInternal <- scaleChecks[[2]]
    
    # =========================================================================
    # COLORS
    # =========================================================================
    
    if (archInternal$clip == TRUE){
        subset <- "pairs_clip"
    } else {
        if (archInternal$clip.noAnchor == TRUE){
            subset <- "pairs_noanchor"
        } else {
            subset <- "pairs"
        }
    }
    
    archColors <- parseColors(data = bedpe,
                                fill = archInternal$fill,
                                object = archesPlot,
                                subset = subset)
    if (length(archColors[[1]]) > 0){
        bedpe$color <- archColors[[1]]
    } else {
        bedpe$color <- rep("#1f4297", nrow(bedpe))
    }
    
    archesPlot <- archColors[[2]]
    
    bedpe$linecolor <- lineColors(linecolor = archInternal$linecolor,
                                    fillcolors = bedpe$color,
                                    data = bedpe,
                                    object = archesPlot,
                                    subset = subset)
    
    # =========================================================================
    # SUBSET DATA
    # =========================================================================

    if (!is.null(archesPlot$chromstart) & !is.null(archesPlot$chromend)) {
        if (archInternal$clip == TRUE) {
            bedpe <- bedpe[which(bedpe[, "chrom1"] == archesPlot$chrom &
                bedpe[, "chrom2"] == archesPlot$chrom &
                bedpe[, "start1"] >= archesPlot$chromstart &
                bedpe[, "end1"] <= archesPlot$chromend &
                bedpe[, "start2"] >= archesPlot$chromstart &
                bedpe[, "end2"] <= archesPlot$chromend), ]
        } else {
            if (archInternal$clip.noAnchor == TRUE){
                bedpe <- bedpe[which(bedpe[, "chrom1"] == archesPlot$chrom &
                                    bedpe[, "chrom2"] == archesPlot$chrom &
                                ((bedpe[, "start1"] >= archesPlot$chromstart &
                                bedpe[, "start1"] <= archesPlot$chromend) |
                                (bedpe[, "end2"] >= archesPlot$chromstart &
                                bedpe[, "end2"] <= archesPlot$chromend))), ]
            } else {
                bedpe <- bedpe[which(bedpe[, "chrom1"] == archesPlot$chrom &
                                        bedpe[, "chrom2"] == archesPlot$chrom),]
                overlappingRanges <- as.data.frame(subsetByOverlaps(ranges = 
                                        IRanges(start = archesPlot$chromstart, 
                                                end = archesPlot$chromend),
                                        x = IRanges(start = bedpe[,"start1"], 
                                                    end = bedpe[,"end2"])))
                bedpe <- bedpe[which(bedpe[,"start1"] %in% 
                                        overlappingRanges$start &
                                        bedpe[,"end2"] %in% 
                                        overlappingRanges$end),]
            }
        }
    } else {
        bedpe <- data.frame(matrix(nrow = 0, ncol = 6))
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(archesPlot$x) | is.null(archesPlot$y)) {
        vp <- viewport(
            height = unit(0.5, "npc"), width = unit(1, "npc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            xscale = archInternal$xscale,
            clip = "on",
            just = "center",
            name = "arches1"
        )

        if (archInternal$draw == TRUE) {
            grid.newpage()
        }
        
    } else {
        
        ## Get viewport name
        currentViewports <- current_viewports()
        vp_name <- paste0(
            "arches",
            length(grep(
                pattern = "arches",
                x = currentViewports
            )) + 1
        )
        
        addViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = archesPlot)

        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            xscale = archInternal$xscale,
            clip = "on",
            just = archInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # HEIGHTS
    # =========================================================================

    if (nrow(bedpe) > 0) {
        if (is.null(archInternal$archHeight)) {
            bedpe$height <- rep(1, nrow(bedpe))
        } else if (is(archInternal$archHeight, "numeric")){
            if (length(archInternal$archHeight) == 1){
                bedpe$height <- rep(archInternal$archHeight, nrow(bedpe))
            } else {
                if (length(archInternal$archHeight) < nrow(bedpe)){
                        stop("`archHeight` vector is shorter than the", 
                            " number of plotted arches.", call. = FALSE)
                } else if (length(archInternal$archHeight) > nrow(bedpe)){
                    warning("`archHeight` vector is longer than the number ",
                            "of plotted arches. `archHeight` vector will",
                            " be truncated.", call. = FALSE)
                }
                bedpe$height <- archInternal$archHeight[seq(1, nrow(bedpe))]
            }
            yscale <- height_yscale(
                        heights = bedpe$height,
                        flip = archInternal$flip,
                        range = archInternal$range
                    )
            vp$yscale <- yscale
        } else {
            if (!archInternal$archHeight %in% colnames(bedpe)){
                stop("Column name for `archHeight` not found in data.",
                    call. = FALSE)
            }
            bedpe$height <- bedpe[, archInternal$archHeight]
            yscale <- height_yscale(
                        heights = bedpe$height,
                        flip = archInternal$flip,
                        range = archInternal$range
                    )
            vp$yscale <- yscale
            
        }

        if (length(bedpe$height) > 0) {
            bedpe$normHeight <- lapply(bedpe$height, normHeights,
                min = min(vp$yscale),
                max = max(vp$yscale)
            )
        }
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================

    backgroundGrob <- rectGrob(
        gp = gpar(fill = archInternal$bg, col = NA),
        name = "background"
    )
    assign("arches_grobs", gTree(vp = vp, children = gList(backgroundGrob)),
        envir = pgEnv
    )

    # =========================================================================
    # GROBS
    # =========================================================================

    if (nrow(bedpe) > 0) {
        if (archInternal$baseline == TRUE) {
            baselineGrob <- segmentsGrob(
                x0 = unit(0, "npc"), y0 = 0,
                x1 = unit(1, "npc"), y1 = 0,
                default.units = "native",
                gp = gpar(
                    col = archInternal$baseline.color,
                    lwd = archInternal$baseline.lwd
                )
            )
            assign("arches_grobs",
                addGrob(
                    gTree = get("arches_grobs", envir = pgEnv),
                    child = baselineGrob
                ),
                envir = pgEnv
            )
        }
        
        invisible(apply(bedpe, 1, drawRibbons,
            style = archInternal$style,
            arch = archInternal$curvature,
            flip = archInternal$flip,
            transp = archInternal$alpha, gp = archInternal$gp
        ))
    } else {
        if (archInternal$txdbChecks == TRUE) {
            warning("Data contains no values.", call. = FALSE)
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (archInternal$draw == TRUE) {
        grid.draw(get("arches_grobs", envir = pgEnv))
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    archesPlot$grobs <- get("arches_grobs", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("arches[", vp$name, "]")
    invisible(archesPlot)
}
