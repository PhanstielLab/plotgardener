#' Plot a Hi-C interaction matrix in a triangular format
#' 
#' @usage bb_plotHicTriangle(
#'     data,
#'     resolution = "auto",
#'     zrange = NULL,
#'     norm = "KR",
#'     matrix = "observed",
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
#'     palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu")),
#'     colorTrans = "linear",
#'     x = NULL,
#'     y = NULL,
#'     width = NULL,
#'     height = NULL,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     draw = TRUE,
#'     params = NULL,
#'     quiet = FALSE
#' )
#'
#' @param data Path to .hic file as a string or a 3-column
#' dataframe of interaction counts in sparse upper triangular format.
#' @param resolution A numeric specifying the width in basepairs of
#' each pixel. For hic files, "auto" will attempt to choose a
#' resolution based on the size of the region. For
#' dataframes, "auto" will attempt to detect the resolution
#' the dataframe contains.
#' @param zrange A numeric vector of length 2 specifying the
#' range of interaction scores to plot, where extreme values
#' will be set to the max or min.
#' @param norm Character value specifying hic data normalization method,
#' if giving .hic file. This value must be found in the .hic file.
#' Default value is \code{norm = "KR"}.
#' @param matrix Character value indicating the type of matrix to output.
#' Default value is \code{matrix = "observed"}. Options are:
#' \itemize{
#' \item{\code{"observed"}: }{Observed counts.}
#' \item{\code{"oe"}: }{Observed/expected counts.}
#' \item{\code{"log2oe"}: }{Log2 transformed observed/expected counts.}
#' }
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[BentoBox]{bb_assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param palette A function describing the color palette to use for
#' representing scale of interaction scores. Default value is
#' \code{palette =  colorRampPalette(brewer.pal(n = 9, "YlGnBu"))}.
#' @param colorTrans A string specifying how to scale Hi-C colors.
#' Options are "linear", "log", "log2", or "log10".
#' Default value is \code{colorTrans = "linear"}.
#' @param x A numeric or unit object specifying triangle Hi-C plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying triangle Hi-C plot y-location.
#' The character value will
#' place the triangle Hi-C plot y relative to the bottom of the most
#' recently plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying the bottom
#' width of the Hi-C plot triangle.
#' @param height A numeric or unit object specifying the height of
#' the Hi-C plot triangle.
#' @param just Justification of triangle Hi-C plot relative to
#' its (x, y) location. If there are two values, the first value specifies
#' horizontal justification and the second value specifies vertical
#' justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if
#' \code{x}, \code{y}, \code{width}, or \code{height} are only given as
#' numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should
#' be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_params} object containing
#' relevant function parameters.
#' @param quiet A logical indicating whether or not to print messages.
#'
#' @return Returns a \code{bb_hicTriangle} object containing relevant
#' genomic region, Hi-C data, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data
#' library(BentoBoxData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create a page
#' bb_pageCreate(width = 4, height = 2.5, default.units = "inches")
#'
#' ## Plot and place triangle Hi-C plot
#' hicPlot <- bb_plotHicTriangle(
#'     data = IMR90_HiC_10kb, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     x = 2, y = 0.5, width = 3, height = 1.5,
#'     just = "top", default.units = "inches"
#' )
#'
#' ## Annotate x-axis genome label
#' bb_annoGenomeLabel(
#'     plot = hicPlot, scale = "Mb", x = 0.5, y = 2.03,
#'     just = c("left", "top")
#' )
#'
#' ## Annotate heatmap legend
#' bb_annoHeatmapLegend(
#'     plot = hicPlot, x = 3.5, y = 0.5,
#'     width = 0.13, height = 1.2,
#'     just = c("right", "top")
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @details
#' A triangle Hi-C plot can be placed on a BentoBox coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' bb_plotHicTriangle(data, chrom,
#'                 chromstart = NULL, chromend = NULL,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated triangle
#' Hi-C plot by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotHicTriangle(data, chrom,
#'                 chromstart = NULL, chromend = NULL)
#' }
#'
#' If \code{height} is \eqn{<} \eqn{0.5 * sqrt(2)}, the top of the triangle
#' will be cropped to the given \code{height}.
#'
#' @seealso \link[BentoBox]{bb_readHic}
#'
#' @export
bb_plotHicTriangle <- function(data, resolution = "auto", zrange = NULL,
                            norm = "KR", matrix = "observed", chrom,
                            chromstart = NULL, chromend = NULL,
                            assembly = "hg38",
                            palette = colorRampPalette(brewer.pal(
                                n = 9, "YlGnBu"
                            )),
                            colorTrans = "linear", x = NULL, y = NULL,
                            width = NULL, height = NULL,
                            just = c("left", "top"),
                            default.units = "inches", draw = TRUE,
                            params = NULL, quiet = FALSE) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for bb_plotTriangleHic
    errorcheck_bb_plotTriangleHic <- function(hic, hic_plot, norm, assembly) {

        ###### hic/norm #####
        bb_hicErrors(hic = hic,
                    norm = norm)

        ## Even though straw technically works without "chr" for hg19,
        ## will not accept for consistency purposes
        if (assembly == "hg19") {
            if (grepl("chr", hic_plot$chrom) == FALSE) {
                stop("'", hic_plot$chrom, "'",
                    "is an invalid input for an hg19 chromsome. ",
                    "Please specify chromosome as",
                    "'chr", hic_plot$chrom, "'.",
                    call. = FALSE
                )
            }
        }
        
        bb_regionErrors(chromstart = hic_plot$chromstart,
                    chromend = hic_plot$chromend)

        ##### zrange #####
        bb_rangeErrors(range = hic_plot$zrange)
        
        ##### height #####
        if (!is.null(hic_plot$height)) {

            ## convert height to inches
            height <- convertHeight(hic_plot$height,
                unitTo = "inches", valueOnly = TRUE
            )
            if (height < 0.05) {
                stop("Height is too small for a valid triangle Hi-C plot.",
                    call. = FALSE
                )
            }
        }
    }

    ## Define a function that subsets data
    subset_data <- function(hic, hic_plot) {
        if (nrow(hic) > 0) {
            hic <- hic[which(hic[, "x"] >= floor(hic_plot$chromstart /
                hic_plot$resolution) *
                hic_plot$resolution &
                hic[, "x"] < hic_plot$chromend &
                hic[, "y"] >= floor(hic_plot$chromstart /
                    hic_plot$resolution) *
                    hic_plot$resolution &
                hic[, "y"] < hic_plot$chromend), ]
        }


        return(hic)
    }

    ## Define a function that manually "clips" squares/triangles along edges
    manual_clip <- function(hic, hic_plot) {
        clipLeft <- hic[which(hic[, "x"] < hic_plot$chromstart), ]
        clipTop <- hic[which((hic[, "y"] + hic_plot$resolution) >
            hic_plot$chromend), ]

        topLeft <- suppressMessages(dplyr::inner_join(clipLeft, clipTop))

        clipLeft <- suppressMessages(dplyr::anti_join(clipLeft, topLeft))
        clipTop <- suppressMessages(dplyr::anti_join(clipTop, topLeft))

        ############# Squares
        squares <- hic[which(hic[, "y"] > hic[, "x"]), ]
        clipLeftsquares <- suppressMessages(dplyr::inner_join(
            squares,
            clipLeft
        ))
        clipTopsquares <- suppressMessages(dplyr::inner_join(
            squares,
            clipTop
        ))
        clippedSquares <- rbind(
            clipLeftsquares, clipTopsquares,
            topLeft
        )

        squares <- suppressMessages(dplyr::anti_join(squares, clippedSquares))


        clipLeftsquares$width <- hic_plot$resolution - (hic_plot$chromstart -
            clipLeftsquares$x)
        clipLeftsquares$x <- rep(hic_plot$chromstart, nrow(clipLeftsquares))

        clipTopsquares$height <- hic_plot$chromend - clipTopsquares$y

        topLeft$width <- hic_plot$resolution - (hic_plot$chromstart -
            topLeft$x)
        topLeft$x <- rep(hic_plot$chromstart, nrow(topLeft))

        topLeft$height <- hic_plot$chromend - topLeft$y


        ############# Triangles
        triangles <- hic[which(hic[, "y"] == hic[, "x"]), ]
        topRight <- suppressMessages(dplyr::inner_join(triangles, clipTop))
        bottomLeft <- suppressMessages(dplyr::inner_join(triangles, clipLeft))
        clippedTriangles <- rbind(topRight, bottomLeft)

        triangles <- suppressMessages(dplyr::anti_join(
            triangles,
            clippedTriangles
        ))

        topRight$height <- hic_plot$chromend - topRight$y
        topRight$width <- topRight$height

        bottomLeft$width <- hic_plot$resolution - (hic_plot$chromstart -
            bottomLeft$x)
        bottomLeft$height <- bottomLeft$width
        bottomLeft$x <- rep(hic_plot$chromstart, nrow(bottomLeft))
        bottomLeft$y <- rep(hic_plot$chromstart, nrow(bottomLeft))


        ## Recombine
        clippedHic <- rbind(
            squares, triangles, clipLeftsquares,
            clipTopsquares, topLeft, topRight, bottomLeft
        )

        return(clippedHic)
    }

    ## Define a function that makes grobs for the triangle hic diagonal
    hic_diagonal <- function(hic) {
        col <- hic["color"]
        x <- utils::type.convert(hic["x"], as.is = TRUE)
        y <- utils::type.convert(hic["y"], as.is = TRUE)
        width <- utils::type.convert(hic["width"], as.is = TRUE)
        height <- utils::type.convert(hic["height"], as.is = TRUE)

        xleft <- x
        xright <- x + width
        ybottom <- y
        ytop <- y + height

        hic_triangle <- polygonGrob(
            x = c(xleft, xleft, xright),
            y = c(ybottom, ytop, ytop),
            gp = gpar(col = NA, fill = col),
            default.units = "native"
        )

        assign("hic_grobs2",
            addGrob(
                gTree = get("hic_grobs2", envir = bbEnv),
                child = hic_triangle
            ),
            envir = bbEnv
        )
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_thicInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_thicInternal"
    )
    ## Justification
    bb_thicInternal$just <- bb_justConversion(just = bb_thicInternal$just)
    
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    hic_plot <- structure(list(
        chrom = bb_thicInternal$chrom,
        chromstart = bb_thicInternal$chromstart,
        chromend = bb_thicInternal$chromend,
        altchrom = bb_thicInternal$chrom,
        altchromstart = bb_thicInternal$chromstart,
        altchromend = bb_thicInternal$chromend,
        assembly = bb_thicInternal$assembly,
        resolution = bb_thicInternal$resolution,
        x = bb_thicInternal$x, y = bb_thicInternal$y,
        width = bb_thicInternal$width,
        height = bb_thicInternal$height,
        just = bb_thicInternal$just,
        color_palette = NULL,
        colorTrans = bb_thicInternal$colorTrans,
        zrange = bb_thicInternal$zrange,
        outsideVP = NULL, grobs = NULL
    ),
    class = "bb_hicTriangle"
    )
    attr(x = hic_plot, which = "plotted") <- bb_thicInternal$draw

    # =========================================================================
    # CHECK PLACEMENT/ARGUMENT ERRORS
    # =========================================================================

    if (is.null(bb_thicInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(bb_thicInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)
    check_placement(object = hic_plot)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    hic_plot$assembly <- parse_bbAssembly(assembly = hic_plot$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    hic_plot <- defaultUnits(
        object = hic_plot,
        default.units = bb_thicInternal$default.units
    )

    # ========================================================================
    # CATCH ERRORS
    # ========================================================================

    errorcheck_bb_plotTriangleHic(
        hic = bb_thicInternal$data,
        hic_plot = hic_plot,
        norm = bb_thicInternal$norm,
        assembly = hic_plot$assembly$Genome
    )

    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================
    
    scaleChecks <- genomicScale(object = hic_plot,
                                objectInternal = bb_thicInternal,
                                plotType = "triangle Hi-C plot")
    hic_plot <- scaleChecks[[1]]
    bb_thicInternal <- scaleChecks[[2]]
    
    # =========================================================================
    # ADJUST RESOLUTION
    # =========================================================================

    if (bb_thicInternal$resolution == "auto") {
        hic_plot <- adjust_resolution(
            hic = bb_thicInternal$data,
            hic_plot = hic_plot
        )
    } else {
        hic_plot <- hic_limit(
            hic = bb_thicInternal$data,
            hic_plot = hic_plot
        )
    }

    # =========================================================================
    # READ IN DATA
    # =========================================================================

    hic <- read_data(
        hic = bb_thicInternal$data, hic_plot = hic_plot,
        norm = bb_thicInternal$norm, assembly = hic_plot$assembly,
        type = bb_thicInternal$matrix,
        quiet = bb_thicInternal$quiet
    )

    # =========================================================================
    # SUBSET DATA
    # =========================================================================

    hic <- subset_data(hic = hic, hic_plot = hic_plot)

    # =========================================================================
    # SET ZRANGE AND SCALE DATA
    # =========================================================================

    hic_plot <- set_zrange(hic = hic, hic_plot = hic_plot)
    hic$counts[hic$counts <= hic_plot$zrange[1]] <- hic_plot$zrange[1]
    hic$counts[hic$counts >= hic_plot$zrange[2]] <- hic_plot$zrange[2]

    # =========================================================================
    # CONVERT NUMBERS TO COLORS
    # =========================================================================

    ## if we don't have an appropriate zrange (even after setting it
    ## based on a null zrange), can't scale to colors
    if (!is.null(hic_plot$zrange) & length(unique(hic_plot$zrange)) == 2) {
        if (grepl("log", bb_thicInternal$colorTrans) == TRUE) {
            logBase <- utils::type.convert(gsub("log", "", 
                                                bb_thicInternal$colorTrans), 
                                as.is = TRUE)
            if (is.na(logBase)) {
                logBase <- exp(1)
            }

            ## Won't scale to log if negative values
            if (any(hic$counts < 0)) {
                stop("Negative values in Hi-C data. Cannot scale colors ",
                "on a log scale. Please ",
                "set `colorTrans = 'linear'`.", call. = FALSE)
            }


            hic$counts <- log(hic$counts, base = logBase)
            hic$color <- bb_mapColors(vector = hic$counts,
                palette = bb_thicInternal$palette,
                range = log(hic_plot$zrange, logBase)
            )
        } else {
            hic$color <- bb_mapColors(hic$counts,
                palette = bb_thicInternal$palette,
                range = hic_plot$zrange
            )
        }

        hic_plot$color_palette <- bb_thicInternal$palette
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_hicTriangle",
        length(grep(
            pattern = "bb_hicTriangle",
            x = currentViewports
        )) + 1
    )

    if (is.null(hic_plot$x) | is.null(hic_plot$y)) {
        inside_vp <- viewport(
            height = unit(1, "npc"), width = unit(0.5, "npc"),
            x = unit(0, "npc"), y = unit(0, "npc"),
            xscale = bb_thicInternal$xscale,
            yscale = bb_thicInternal$xscale,
            just = c("left", "bottom"),
            name = paste0(vp_name, "_inside"),
            angle = -45
        )

        outside_vp <- viewport(
            height = unit(0.75, "snpc"),
            width = unit(1.5, "snpc"),
            x = unit(0.125, "npc"),
            y = unit(0.25, "npc"),
            xscale = bb_thicInternal$xscale,
            clip = "on",
            just = c("left", "bottom"),
            name = paste0(vp_name, "_outside")
        )


        if (bb_thicInternal$draw == TRUE) {
            vp_name <- "bb_hicTriangle1"
            inside_vp$name <- "bb_hicTriangle1_inside"
            outside_vp$name <- "bb_hicTriangle1_outside"
            grid.newpage()
        }
    } else {

        ## Get sides of viewport based on input width
        vp_side <- (convertWidth(hic_plot$width,
            unitTo = get("page_units", envir = bbEnv),
            valueOnly = TRUE
        )) / sqrt(2)

        ## Convert coordinates into same units as page for outside vp
        page_coords <- convert_page(object = hic_plot)

        ## Get bottom left point of triangle (hence bottom left of actual
        ## viewport) based on just
        bottom_coords <- vp_bottomLeft(viewport(
            x = page_coords$x,
            y = page_coords$y,
            width = page_coords$width,
            height = page_coords$height,
            just = hic_plot$just
        ))

        inside_vp <- viewport(
            height = unit(
                vp_side,
                get("page_units", envir = bbEnv)
            ),
            width = unit(
                vp_side,
                get("page_units", envir = bbEnv)
            ),
            x = unit(0, "npc"),
            y = unit(0, "npc"),
            xscale = bb_thicInternal$xscale,
            yscale = bb_thicInternal$xscale,
            just = c("left", "bottom"),
            name = paste0(vp_name, "_inside"),
            angle = -45
        )

        outside_vp <- viewport(
            height = page_coords$height,
            width = page_coords$width,
            x = unit(
                bottom_coords[[1]],
                get("page_units", envir = bbEnv)
            ),
            y = unit(
                bottom_coords[[2]],
                get("page_units", envir = bbEnv)
            ),
            xscale = bb_thicInternal$xscale,
            clip = "on",
            just = c("left", "bottom"),
            name = paste0(vp_name, "_outside")
        )

        add_bbViewport(paste0(vp_name, "_outside"))
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS
    # =========================================================================

    hic_plot$outsideVP <- outside_vp
    assign("hic_grobs2", gTree(vp = inside_vp), envir = bbEnv)

    # =========================================================================
    # MAKE GROBS
    # =========================================================================

    if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)) {
        if (nrow(hic) > 0) {
            hic$width <- hic_plot$resolution
            hic$height <- hic_plot$resolution

            ## Manually "clip" the grobs that fall out of the desired chmstart
            ## to chromend region
            hic <- manual_clip(hic = hic, hic_plot = hic_plot)
            hic <- hic[order(utils::type.convert(rownames(hic), 
                                            as.is = TRUE)), ]

            ## Separate into squares for upper region and triangle
            ## shapes for the diagonal
            squares <- hic[which(hic[, "y"] > hic[, "x"]), ]
            triangles <- hic[which(hic[, "y"] == hic[, "x"]), ]

            if (nrow(squares) > 0) {

                ## Make square grobs and add to grob gTree
                hic_squares <- rectGrob(
                    x = squares$x,
                    y = squares$y,
                    just = c("left", "bottom"),
                    width = squares$width,
                    height = squares$height,
                    gp = gpar(col = NA, fill = squares$color),
                    default.units = "native"
                )
                assign("hic_grobs2",
                    addGrob(
                        gTree = get("hic_grobs2", envir = bbEnv),
                        child = hic_squares
                    ),
                    envir = bbEnv
                )
            }

            if (nrow(triangles) > 0) {
                ## Make triangle grobs and add to grob gTree
                invisible(apply(triangles, 1, hic_diagonal))
            }

            if (nrow(squares) == 0 & nrow(triangles) == 0) {
                if (bb_thicInternal$txdbChecks == TRUE) {
                    if (!is.na(hic_plot$resolution)){
                        warning("No data found in region. Suggestions: ",
                                "check that ",
                                "chromosome names match genome assembly; ",
                                "check region.", call. = FALSE)
                    }
                }
            }
        } else {
            if (!is.na(hic_plot$resolution)){
                warning("No data found in region. Suggestions: check that ",
                        "chromosome names match genome assembly; ",
                        "check region.", call. = FALSE)
            }
        }
    }


    # =========================================================================
    # IF DRAW == TRUE, DRAW GROBS
    # =========================================================================

    if (bb_thicInternal$draw == TRUE) {
        pushViewport(outside_vp)
        grid.draw(get("hic_grobs2", envir = bbEnv))
        upViewport()
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    hic_plot$grobs <- get("hic_grobs2", envir = bbEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_hicTriangle[", vp_name, "]")
    invisible(hic_plot)
}
