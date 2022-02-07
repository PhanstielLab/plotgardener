#' Plot a Hi-C interaction matrix in a triangular format
#' 
#' @usage plotHicTriangle(
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
#'     flip = FALSE,
#'     bg = NA,
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
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param palette A function describing the color palette to use for
#' representing scale of interaction scores. Default value is
#' \code{palette =  colorRampPalette(brewer.pal(n = 9, "YlGnBu"))}.
#' @param colorTrans A string specifying how to scale Hi-C colors.
#' Options are "linear", "log", "log2", or "log10".
#' Default value is \code{colorTrans = "linear"}.
#' @param flip A logical indicating whether to flip the orientation of
#' the Hi-C matrix over the x-axis. Default value is \code{flip = FALSE}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param x A numeric or unit object specifying triangle Hi-C plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying triangle Hi-C plot y-location.
#' The character value will
#' place the triangle Hi-C plot y relative to the bottom of the most
#' recently plotted plot according to the units of the plotgardener page.
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
#' @param params An optional \link[plotgardener]{pgParams} object containing
#' relevant function parameters.
#' @param quiet A logical indicating whether or not to print messages.
#'
#' @return Returns a \code{hicTriangle} object containing relevant
#' genomic region, Hi-C data, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create a page
#' pageCreate(width = 4, height = 2.5, default.units = "inches")
#'
#' ## Plot and place triangle Hi-C plot
#' hicPlot <- plotHicTriangle(
#'     data = IMR90_HiC_10kb, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19", bg = "black",
#'     x = 2, y = 0.5, width = 3, height = 1.5,
#'     just = "top", default.units = "inches"
#' )
#'
#' ## Annotate x-axis genome label
#' annoGenomeLabel(
#'     plot = hicPlot, scale = "Mb", x = 0.5, y = 2.03,
#'     just = c("left", "top")
#' )
#'
#' ## Annotate heatmap legend
#' annoHeatmapLegend(
#'     plot = hicPlot, x = 3.5, y = 0.5,
#'     width = 0.13, height = 1.2,
#'     just = c("right", "top")
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' A triangle Hi-C plot can be placed on a plotgardener coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' plotHicTriangle(data, chrom,
#'                 chromstart = NULL, chromend = NULL,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated triangle
#' Hi-C plot by ignoring plot placement parameters:
#' \preformatted{
#' plotHicTriangle(data, chrom,
#'                 chromstart = NULL, chromend = NULL)
#' }
#'
#' If \code{height} is \eqn{<} \eqn{0.5 * width}, the top of the triangle
#' will be cropped to the given \code{height}.
#'
#' @seealso \link[plotgardener]{readHic}
#'
#' @export
plotHicTriangle <- function(data, resolution = "auto", zrange = NULL,
                            norm = "KR", matrix = "observed", chrom,
                            chromstart = NULL, chromend = NULL,
                            assembly = "hg38",
                            palette = colorRampPalette(brewer.pal(
                                n = 9, "YlGnBu"
                            )),
                            colorTrans = "linear", flip = FALSE, 
                            bg = NA,
                            x = NULL, y = NULL,
                            width = NULL, height = NULL,
                            just = c("left", "top"),
                            default.units = "inches", draw = TRUE,
                            params = NULL, quiet = FALSE) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for plotTriangleHic
    errorcheck_plotTriangleHic <- function(hic, hicPlot, norm) {

        ###### hic/norm #####
        hicErrors(hic = hic,
                    norm = norm)

        
        regionErrors(chromstart = hicPlot$chromstart,
                    chromend = hicPlot$chromend)

        ##### zrange #####
        rangeErrors(range = hicPlot$zrange)
        
        ##### height #####
        if (!is.null(hicPlot$height)) {

            ## convert height to inches
            height <- convertHeight(hicPlot$height,
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
    subset_data <- function(hic, hicPlot) {
        if (nrow(hic) > 0) {
            hic <- hic[which(hic[, "x"] >= floor(hicPlot$chromstart /
                hicPlot$resolution) *
                hicPlot$resolution &
                hic[, "x"] < hicPlot$chromend &
                hic[, "y"] >= floor(hicPlot$chromstart /
                    hicPlot$resolution) *
                    hicPlot$resolution &
                hic[, "y"] < hicPlot$chromend), ]
        }


        return(hic)
    }

    ## Define a function that manually "clips" pixels along edges
    manual_clip <- function(hic, hicPlot, flip) {
        
        if (flip == FALSE){
            clipLeft <- hic[which(hic[, "x"] < hicPlot$chromstart), ]
            clipTop <- hic[which((hic[, "y"] + hicPlot$resolution) >
                hicPlot$chromend), ]

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

            squares <- suppressMessages(dplyr::anti_join(squares, 
                                                        clippedSquares))


            clipLeftsquares$width <- hicPlot$resolution - (hicPlot$chromstart -
                clipLeftsquares$x)
            clipLeftsquares$x <- rep(hicPlot$chromstart, nrow(clipLeftsquares))

            clipTopsquares$height <- hicPlot$chromend - clipTopsquares$y

            topLeft$width <- hicPlot$resolution - (hicPlot$chromstart -
                topLeft$x)
            topLeft$x <- rep(hicPlot$chromstart, nrow(topLeft))

            topLeft$height <- hicPlot$chromend - topLeft$y


            ############# Triangles
            triangles <- hic[which(hic[, "y"] == hic[, "x"]), ]
            topRight <- suppressMessages(dplyr::inner_join(triangles, clipTop))
            bottomLeft <- suppressMessages(dplyr::inner_join(triangles, 
                                                                    clipLeft))
            clippedTriangles <- rbind(topRight, bottomLeft)

            triangles <- suppressMessages(dplyr::anti_join(
                triangles,
                clippedTriangles
            ))

            topRight$height <- hicPlot$chromend - topRight$y
            topRight$width <- topRight$height

            bottomLeft$width <- hicPlot$resolution - (hicPlot$chromstart -
                bottomLeft$x)
            bottomLeft$height <- bottomLeft$width
            bottomLeft$x <- rep(hicPlot$chromstart, nrow(bottomLeft))
            bottomLeft$y <- rep(hicPlot$chromstart, nrow(bottomLeft))
            
            clippedHic <- rbind(
                    squares, triangles, clipLeftsquares,
                    clipTopsquares, topLeft, topRight, bottomLeft
                )
            
        } else {
            
            clipBottom <- hic[which(hic[, "y"] < hicPlot$chromstart), ]
            clipRight <- hic[which((hic[, "x"] + hicPlot$resolution) >
                                hicPlot$chromend), ]
            
            bottomRight <- suppressMessages(dplyr::inner_join(clipBottom, 
                                                                clipRight))
            
            clipBottom <- suppressMessages(dplyr::anti_join(clipBottom, 
                                                            bottomRight))
            clipRight <- suppressMessages(dplyr::anti_join(clipRight, 
                                                            bottomRight))
            
            ############# Squares
            squares <- hic[which(hic[, "x"] > hic[, "y"]), ]
            clipBottomsquares <- suppressMessages(dplyr::inner_join(
                squares,
                clipBottom
            ))
            clipRightsquares <- suppressMessages(dplyr::inner_join(
                squares,
                clipRight
            ))
            clippedSquares <- rbind(
                clipBottomsquares, clipRightsquares,
                bottomRight
            )
            
            squares <- suppressMessages(dplyr::anti_join(squares, 
                                                        clippedSquares))
            
            clipRightsquares$width <- hicPlot$chromend - clipRightsquares$x
            clipBottomsquares$height <- hicPlot$resolution - 
                (hicPlot$chromstart - clipBottomsquares$y)
            clipBottomsquares$y <- rep(hicPlot$chromstart, 
                                        nrow(clipBottomsquares))
            
            
            bottomRight$height <- hicPlot$resolution - (hicPlot$chromstart -
                                        bottomRight$y)
            bottomRight$y <- rep(hicPlot$chromstart, nrow(bottomRight))
            bottomRight$width <- hicPlot$chromend - bottomRight$x
            
            
            ############# Triangles
            triangles <- hic[which(hic[, "y"] == hic[, "x"]), ]
            topRight <- suppressMessages(dplyr::inner_join(triangles, 
                                                            clipRight))
            bottomLeft <- suppressMessages(dplyr::inner_join(triangles, 
                                                            clipBottom))
            
            
            clippedTriangles <- rbind(topRight, bottomLeft)
            
            triangles <- suppressMessages(dplyr::anti_join(
                triangles,
                clippedTriangles
            ))
            
            topRight$width <- hicPlot$chromend - topRight$x
            topRight$height <- topRight$width
            
            bottomLeft$height <- hicPlot$resolution - (hicPlot$chromstart -
                                        bottomLeft$y)
            bottomLeft$width <- bottomLeft$height
            bottomLeft$x <- rep(hicPlot$chromstart, nrow(bottomLeft))
            bottomLeft$y <- rep(hicPlot$chromstart, nrow(bottomLeft))
            
            
            clippedHic <- rbind(
                squares, triangles, clipRightsquares,
                clipBottomsquares, bottomRight, topRight, bottomLeft
            )
            
        }

        return(clippedHic)
    }

    ## Define a function that makes grobs for the triangle hic diagonal
    hic_diagonal <- function(hic, flip) {
        col <- hic["color"]
        x <- utils::type.convert(hic["x"], as.is = TRUE)
        y <- utils::type.convert(hic["y"], as.is = TRUE)
        width <- utils::type.convert(hic["width"], as.is = TRUE)
        height <- utils::type.convert(hic["height"], as.is = TRUE)
        
        xleft <- x
        xright <- x + width
        ybottom <- y
        ytop <- y + height
        
        if (flip == FALSE){

            hic_triangle <- polygonGrob(
                x = c(xleft, xleft, xright),
                y = c(ybottom, ytop, ytop),
                gp = gpar(col = NA, fill = col),
                default.units = "native"
            )
            
        } else {
            
            hic_triangle <- polygonGrob(
                x = c(xleft, xright, xright),
                y = c(ybottom, ybottom, ytop),
                gp = gpar(col = NA, fill = col),
                default.units = "native"
            )
        }

        assign("hic_grobs2",
            addGrob(
                gTree = get("hic_grobs2", envir = pgEnv),
                child = hic_triangle
            ),
            envir = pgEnv
        )
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    thicInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "thicInternal"
    )
    ## Justification
    thicInternal$just <- justConversion(just = thicInternal$just)
    
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    hicPlot <- structure(list(
        chrom = thicInternal$chrom,
        chromstart = thicInternal$chromstart,
        chromend = thicInternal$chromend,
        altchrom = thicInternal$chrom,
        altchromstart = thicInternal$chromstart,
        altchromend = thicInternal$chromend,
        assembly = thicInternal$assembly,
        resolution = thicInternal$resolution,
        x = thicInternal$x, y = thicInternal$y,
        width = thicInternal$width,
        height = thicInternal$height,
        just = thicInternal$just,
        color_palette = NULL,
        colorTrans = thicInternal$colorTrans,
        zrange = thicInternal$zrange,
        outsideVP = NULL, grobs = NULL
    ),
    class = "hicTriangle"
    )
    attr(x = hicPlot, which = "plotted") <- thicInternal$draw

    # =========================================================================
    # CHECK PLACEMENT/ARGUMENT ERRORS
    # =========================================================================

    if (is.null(thicInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(thicInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)
    check_placement(object = hicPlot)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    hicPlot$assembly <- parseAssembly(assembly = hicPlot$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    hicPlot <- defaultUnits(
        object = hicPlot,
        default.units = thicInternal$default.units
    )

    # ========================================================================
    # CATCH ERRORS
    # ========================================================================

    errorcheck_plotTriangleHic(
        hic = thicInternal$data,
        hicPlot = hicPlot,
        norm = thicInternal$norm
    )

    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================
    
    scaleChecks <- genomicScale(object = hicPlot,
                                objectInternal = thicInternal,
                                plotType = "triangle Hi-C plot")
    hicPlot <- scaleChecks[[1]]
    thicInternal <- scaleChecks[[2]]
    
    # =========================================================================
    # ADJUST RESOLUTION
    # =========================================================================

    if (thicInternal$resolution == "auto") {
        hicPlot <- adjust_resolution(
            hic = thicInternal$data,
            hicPlot = hicPlot
        )
    } else {
        hicPlot <- hic_limit(
            hic = thicInternal$data,
            hicPlot = hicPlot
        )
    }

    # =========================================================================
    # READ IN DATA
    # =========================================================================

    hic <- read_data(
        hic = thicInternal$data, hicPlot = hicPlot,
        norm = thicInternal$norm, assembly = hicPlot$assembly,
        type = thicInternal$matrix,
        quiet = thicInternal$quiet
    )

    # =========================================================================
    # SUBSET DATA
    # =========================================================================

    hic <- subset_data(hic = hic, hicPlot = hicPlot)
    
    # =========================================================================
    # GET LOWER TRIANGULAR FOR FLIP
    # =========================================================================
    
    if (thicInternal$flip == TRUE){
        hic <- hic[, c("y", "x", "counts")]
        colnames(hic) <- c("x", "y", "counts")
    }

    # =========================================================================
    # SET ZRANGE AND SCALE DATA
    # =========================================================================

    hicPlot <- set_zrange(hic = hic, hicPlot = hicPlot)
    hic$counts[hic$counts <= hicPlot$zrange[1]] <- hicPlot$zrange[1]
    hic$counts[hic$counts >= hicPlot$zrange[2]] <- hicPlot$zrange[2]

    # =========================================================================
    # CONVERT NUMBERS TO COLORS
    # =========================================================================

    ## if we don't have an appropriate zrange (even after setting it
    ## based on a null zrange), can't scale to colors
    if (!is.null(hicPlot$zrange) & length(unique(hicPlot$zrange)) == 2) {
        if (grepl("log", thicInternal$colorTrans) == TRUE) {
            logBase <- utils::type.convert(gsub("log", "", 
                                                thicInternal$colorTrans), 
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
            hic$color <- mapColors(vector = hic$counts,
                palette = thicInternal$palette,
                range = log(hicPlot$zrange, logBase)
            )
        } else {
            hic$color <- mapColors(hic$counts,
                palette = thicInternal$palette,
                range = hicPlot$zrange
            )
        }

        hicPlot$color_palette <- thicInternal$palette
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "hicTriangle",
        length(grep(
            pattern = "hicTriangle",
            x = currentViewports
        )) + 1
    )
    
    ## Inner viewport y-coordinate
    if (thicInternal$flip == TRUE){
        y <- unit(1, "npc")
        
    } else {
        y <- unit(0, "npc")
    }

    if (is.null(hicPlot$x) | is.null(hicPlot$y)) {
        inside_vp <- viewport(
            height = unit(1, "npc"), width = unit(0.5, "npc"),
            x = unit(0, "npc"), 
            y = y,
            xscale = thicInternal$xscale,
            yscale = thicInternal$xscale,
            just = c("left", "bottom"),
            name = paste0(vp_name, "_inside"),
            angle = -45
        )

        outside_vp <- viewport(
            height = unit(0.75, "snpc"),
            width = unit(1.5, "snpc"),
            x = unit(0.25, "npc"),
            y = unit(0.125, "npc"),
            xscale = thicInternal$xscale,
            clip = "on",
            just = c("left", "bottom"),
            name = paste0(vp_name, "_outside")
        )

        if (thicInternal$draw == TRUE) {
            vp_name <- "hicTriangle1"
            inside_vp$name <- "hicTriangle1_inside"
            outside_vp$name <- "hicTriangle1_outside"
            grid.newpage()
        }
    } else {

        ## Get sides of viewport based on input width
        vp_side <- (convertWidth(hicPlot$width,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )) / sqrt(2)

        ## Convert coordinates into same units as page for outside vp
        page_coords <- convert_page(object = hicPlot)

        ## Get bottom left point of triangle (hence bottom left of actual
        ## viewport) based on just
        bottom_coords <- vp_bottomLeft(viewport(
            x = page_coords$x,
            y = page_coords$y,
            width = page_coords$width,
            height = page_coords$height,
            just = hicPlot$just
        ))

        inside_vp <- viewport(
            height = unit(
                vp_side,
                get("page_units", envir = pgEnv)
            ),
            width = unit(
                vp_side,
                get("page_units", envir = pgEnv)
            ),
            x = unit(0, "npc"),
            y = y,
            xscale = thicInternal$xscale,
            yscale = thicInternal$xscale,
            just = c("left", "bottom"),
            name = paste0(vp_name, "_inside"),
            angle = -45
        )

        outside_vp <- viewport(
            height = page_coords$height,
            width = page_coords$width,
            x = unit(
                bottom_coords[[1]],
                get("page_units", envir = pgEnv)
            ),
            y = unit(
                bottom_coords[[2]],
                get("page_units", envir = pgEnv)
            ),
            xscale = thicInternal$xscale,
            clip = "on",
            just = c("left", "bottom"),
            name = paste0(vp_name, "_outside")
        )

        addViewport(paste0(vp_name, "_outside"))
    }

    
    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================
    hicPlot$outsideVP <- outside_vp
    backgroundGrob <- rectGrob(gp = gpar(
        fill = thicInternal$bg,
        col = NA
    ), name = "background")
    assign("hic_grobs2", gTree(vp = inside_vp, children = gList(backgroundGrob)),
        envir = pgEnv
    )
    
    # =========================================================================
    # MAKE GROBS
    # =========================================================================

    if (!is.null(hicPlot$chromstart) & !is.null(hicPlot$chromend)) {
        if (nrow(hic) > 0) {
            
            hic$width <- hicPlot$resolution
            hic$height <- hicPlot$resolution

            ## Manually "clip" the grobs that fall out of the desired chmstart
            ## to chromend region
            hic <- manual_clip(hic = hic, hicPlot = hicPlot, 
                                flip = thicInternal$flip)
            hic <- hic[order(utils::type.convert(rownames(hic),
                                            as.is = TRUE)), ]

            ## Separate into squares for upper region and triangle
            ## shapes for the diagonal
            
            if (thicInternal$flip == TRUE){
                squares <- hic[which(hic[, "x"] > hic[, "y"]), ] 
            } else {
                squares <- hic[which(hic[, "y"] > hic[, "x"]), ]  
            }
            
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
                        gTree = get("hic_grobs2", envir = pgEnv),
                        child = hic_squares
                    ),
                    envir = pgEnv
                )
            }

            if (nrow(triangles) > 0) {
                ## Make triangle grobs and add to grob gTree
                invisible(apply(triangles, 1, hic_diagonal, 
                                flip = thicInternal$flip))
            }

            if (nrow(squares) == 0 & nrow(triangles) == 0) {
                if (thicInternal$txdbChecks == TRUE) {
                    if (!is.na(hicPlot$resolution)){
                        warning("No data found in region.", call. = FALSE)
                    }
                }
            }
        } else {
            if (!is.na(hicPlot$resolution)){
                warning("No data found in region.", call. = FALSE)
            }
        }
    }


    # =========================================================================
    # IF DRAW == TRUE, DRAW GROBS
    # =========================================================================

    if (thicInternal$draw == TRUE) {
        pushViewport(outside_vp)
        
        grid.draw(get("hic_grobs2", envir = pgEnv))
        upViewport()
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    hicPlot$grobs <- get("hic_grobs2", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("hicTriangle[", vp_name, "]")
    invisible(hicPlot)
}
