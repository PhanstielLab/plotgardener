#' Plot a Hi-C interaction matrix in a square format
#' 
#' @usage plotHicSquare(
#'     data,
#'     resolution = "auto",
#'     zrange = NULL,
#'     norm = "KR",
#'     matrix = "observed",
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     altchrom = NULL,
#'     altchromstart = NULL,
#'     altchromend = NULL,
#'     assembly = "hg38",
#'     palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu")),
#'     colorTrans = "linear",
#'     half = "both",
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
#' @param data Path to .hic file as a string or a 3-column dataframe of
#' interaction counts in sparse upper triangular format.
#' @param resolution A numeric specifying the width in
#' basepairs of each pixel. For hic files, "auto" will attempt
#' to choose a resolution based on the size of the region. For
#' dataframes, "auto" will attempt to detect the resolution the
#' dataframe contains.
#' @param zrange A numeric vector of length 2 specifying the range
#' of interaction scores to plot, where extreme values will be set
#' to the max or min.
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
#' @param altchrom Alternate chromosome for off-diagonal plotting
#' or interchromosomal plotting, as a string.
#' @param altchromstart Alternate chromosome integer start position
#' for off-diagonal plotting or interchromosomal plotting.
#' @param altchromend Alternate chromosome integer end position
#' for off-diagonal plotting or interchromosomal plotting.
#' @param assembly Default genome assembly as a string or a
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param palette A function describing the color palette to use for
#' representing scale of interaction scores.
#' Default value is
#' \code{palette =  colorRampPalette(brewer.pal(n = 9, "YlGnBu"))}.
#' @param colorTrans A string specifying how to scale Hi-C colors.
#' Options are "linear", "log", "log2", or "log10".
#' Default value is \code{colorTrans = "linear"}.
#' @param half A character value indicating which diagonal regions to plot.
#' For intrachromosomal plotting, options are \code{"both"}, \code{"top"},
#' or \code{"bottom"}. For off-diagonal or interchromosomal plotting,
#' options are \code{"top"} or \code{"bottom"}.
#' Default value is \code{half = "both"}.
#' \itemize{
#' \item{\code{"both"}: }{Both diagonal halves.}
#' \item{\code{"top"}: }{Half above the diagonal.}
#' \item{\code{"bottom"}: }{Half below the diagonal.}
#' }
#' @param x A numeric or unit object specifying square Hi-C plot x-location.
#' @param y A numeric, unit object, or character containing a "b" combined
#' with a numeric value specifying square Hi-C plot y-location.
#' The character value will
#' place the square Hi-C plot y relative to the bottom of the most recently
#' plotted plot according to the units of the plotgardener page.
#' @param width A numeric or unit object specifying square Hi-C plot width.
#' @param height A numeric or unit object specifying square Hi-C plot height.
#' @param just Justification of square Hi-C plot relative to
#' its (x, y) location. If there are two values, the first value specifies
#' horizontal justification and the second value specifies vertical
#' justification.
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
#' @param quiet A logical indicating whether or not to print messages.
#'
#' @return Returns a \code{hicSquare} object containing relevant
#' genomic region, Hi-C data, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create a page
#' pageCreate(width = 3, height = 3, default.units = "inches")
#'
#' ## Plot and place Hi-C plot
#' hicPlot <- plotHicSquare(
#'     data = IMR90_HiC_10kb, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     x = 0.5, y = 0.5, width = 2, height = 2,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate heatmap legend
#' annoHeatmapLegend(
#'     plot = hicPlot, x = 2.6, y = 0.5,
#'     width = 0.12, height = 1.2,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate x-axis and y-axis genome labels
#' annoGenomeLabel(
#'     plot = hicPlot, scale = "Mb", axis = "x",
#'     x = 0.5, y = 2.53, just = c("left", "top")
#' )
#' annoGenomeLabel(
#'     plot = hicPlot, scale = "Mb", axis = "y",
#'     x = 0.47, y = 0.5, just = c("right", "top")
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' A square Hi-C plot can be placed on a plotgardener coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' plotHicSquare(data, chrom,
#'                 chromstart = NULL, chromend = NULL,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can be used to quickly plot an unannotated
#' square Hi-C plot by ignoring plot placement parameters:
#' \preformatted{
#' plotHicSquare(data, chrom,
#'                 chromstart = NULL, chromend = NULL)
#' }
#'
#' @seealso \link[plotgardener]{readHic}
#'
#' @export
plotHicSquare <- function(data, resolution = "auto", zrange = NULL,
                            norm = "KR", matrix = "observed", chrom,
                            chromstart = NULL, chromend = NULL,
                            altchrom = NULL,
                            altchromstart = NULL, altchromend = NULL,
                            assembly = "hg38",
                            palette = colorRampPalette(brewer.pal(
                                n = 9, "YlGnBu"
                            )),
                            colorTrans = "linear",
                            half = "both", x = NULL, y = NULL, width = NULL,
                            height = NULL, just = c("left", "top"),
                            default.units = "inches",
                            draw = TRUE, params = NULL, quiet = FALSE) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for plotHicSquare
    errorcheck_plotHicSquare <- function(hic, hicPlot, norm) {
        
        ###### hic/norm 
        hicErrors(hic = hic, norm = norm)

        ###### chrom/chromstart/chromend/
        ###### altchrom/altchromstart/altchromend

        ## Even though straw technically works without "chr" for hg19,
        ## will not accept for consistency purposes
        if (hicPlot$assembly$Genome == "hg19") {
            if (grepl("chr", hicPlot$chrom) == FALSE) {
                stop("'", hicPlot$chrom, "'",
                    "is an invalid input for an hg19 chromsome. ",
                    "Please specify chromosome as",
                    "'chr", hicPlot$chrom, "'.",
                    call. = FALSE
                )
            }
        }
        
        regionErrors(chromstart = hicPlot$chromstart,
                    chromend = hicPlot$chromend)
        
        if (!is.null(hicPlot$altchrom)) {

            ## Can't specify altchrom without a chrom
            if (is.null(hicPlot$chrom)) {
                stop("Specified \'altchrom\', but did not give \'chrom\'.",
                    call. = FALSE
                )
            }
            ## Even though straw technically works without "chr" for hg19,
            ## will not accept for consistency purposes
            if (hicPlot$assembly$Genome == "hg19") {
                if (grepl("chr", hicPlot$altchrom) == FALSE) {
                    stop("'", hicPlot$altchrom, "'",
                        "is an invalid input for an hg19 chromsome. ",
                        "Please specify chromosome as",
                        "'chr", hicPlot$altchrom, "'.",
                        call. = FALSE
                    )
                }
            }

            regionErrors(chromstart = hicPlot$altchromstart,
                        chromend = hicPlot$altchromend)


            if (!is.null(hicPlot$chromstart) & !is.null(hicPlot$chromend) &
                !is.null(hicPlot$altchromstart) &
                !is.null(hicPlot$altchromend)) {

                ## Check to see if region is square
                if ((hicPlot$chromend - hicPlot$chromstart) !=
                    (hicPlot$altchromend - hicPlot$altchromstart)) {
                    warning("Trying to plot non-square region.", call. = FALSE)
                }
            }
        }

        ###### zrange #####
        rangeErrors(range = hicPlot$zrange)

        ###### half #####

        if (is.null(hicPlot$altchrom)) {
            if (!(hicPlot$half %in% c("both", "top", "bottom"))) {
                stop("Invalid \'half\'.  Options are \'both\', ",
                    "\top\', or \'bottom\'.",
                    call. = FALSE
                )
            }
        } else {
            if (!hicPlot$half %in% c("top", "bottom")) {
                stop("Invalid \'half\' for off-diagonal and ",
                    "interchromosomal plotting.  Options are \'top\' ",
                    "or \'bottom\'.", call. = FALSE)
            }
        }
        
    }

    ## Define a function that checks for and gets whole chromosome
    ## starts/ends based on an assembly TxDb
    get_wholeChrom <- function(chrom, assembly) {
        chromstart <- NULL
        chromend <- NULL

        if (is(assembly$TxDb, "TxDb")) {
            txdbChecks <- TRUE
        } else {
            
            if (!requireNamespace(assembly$TxDb, quietly = TRUE)){
                txdbChecks <- FALSE
                warning("`", assembly$TxDb, "` not available. Please install",
                        " to plot full chromosome Hi-C map.", call. = FALSE)
            } else {
                txdbChecks <- TRUE
            }
            
        }
        if (txdbChecks == TRUE) {
            if (is(assembly$TxDb, "TxDb")) {
                tx_db <- assembly$TxDb
            } else {
                tx_db <- eval(parse(text = paste0(as.name(assembly$TxDb),
                                                "::",
                                                as.name(assembly$TxDb))))
            }

            assembly_data <- GenomeInfoDb::seqlengths(tx_db)

            if (!chrom %in% names(assembly_data)) {
                warning("Chromosome",
                    "'", chrom, "'",
                    "not found in", "`", assembly$TxDb$packageName, "`",
                    "and data for entire chromosome cannot be plotted.",
                    call. = FALSE
                )
            } else {
                chromstart <- 1
                chromend <- assembly_data[[chrom]]
            }
        }


        return(list(chromstart, chromend))
    }

    ## Define a function that subsets data
    subset_data <- function(hic, hicPlot) {
        if (nrow(hic) > 0) {
            if (is.null(hicPlot$altchrom)) {
                hic <- hic[which(hic[, "x"] >= hicPlot$chromstart -
                    hicPlot$resolution &
                    hic[, "x"] <= hicPlot$chromend + hicPlot$resolution &
                    hic[, "y"] >= hicPlot$chromstart - hicPlot$resolution &
                    hic[, "y"] <= hicPlot$chromend + hicPlot$resolution), ]
            } else {
                if (hicPlot$chrom != hicPlot$altchrom) {
                    hic <- hic[which(hic[, "x"] >= hicPlot$chromstart -
                        hicPlot$resolution &
                        hic[, "x"] <= hicPlot$chromend + hicPlot$resolution &
                        hic[, "y"] >= hicPlot$altchromstart -
                            hicPlot$resolution &
                        hic[, "y"] <= hicPlot$altchromend +
                            hicPlot$resolution), ]
                }
            }
        }

        return(hic)
    }

    ## Define a function that sets viewport xscale and yscale
    vp_scale <- function(hicPlot) {
        if (!is.null(hicPlot$chromstart) & !is.null(hicPlot$chromend)) {
            if (is.null(hicPlot$altchrom)) {
                xscale <- c(hicPlot$chromstart, hicPlot$chromend)
                yscale <- xscale
            } else {
                if (hicPlot$half == "bottom") {
                    xscale <- c(hicPlot$chromstart, hicPlot$chromend)
                    yscale <- c(hicPlot$altchromstart, hicPlot$altchromend)
                }

                if (hicPlot$half == "top") {
                    xscale <- c(hicPlot$altchromstart, hicPlot$altchromend)
                    yscale <- c(hicPlot$chromstart, hicPlot$chromend)
                }
            }
        } else {
            xscale <- c(0, 1)
            yscale <- c(0, 1)
        }

        return(list(xscale, yscale))
    }

    ## Define a function that subsets the hic dataframe into which will
    ## be squares and which will be triangles
    hic_shapes <- function(hic, hicPlot, half) {
        if (!is.null(hicPlot$altchrom)) {

            ## all squares
            squares <- hic
            triangles <- NULL
        } else {
            if (half == "both") {

                ## all squares
                squares <- hic
                triangles <- NULL
            } else if (half == "top") {

                ## squares for top half
                ## triangles for diagonal

                squares <- hic[which(hic[, "y"] > hic[, "x"]), ]
                triangles <- hic[which(hic[, "y"] == hic[, "x"]), ]
            } else if (half == "bottom") {

                ## squares for bottom half
                ## triangles for diagonal

                squares <- hic[which(hic[, "y"] < hic[, "x"]), ]
                triangles <- hic[which(hic[, "y"] == hic[, "x"]), ]
            }
        }

        return(list(squares, triangles))
    }

    ## Define a function that makes grobs for the square hic diagonal
    hic_diagonal <- function(hic, hicPlot, half) {
        col <- hic["color"]
        x <- utils::type.convert(hic["x"], as.is = TRUE)
        y <- utils::type.convert(hic["y"], as.is = TRUE)

        xleft <- x
        xright <- x + hicPlot$resolution
        ybottom <- y
        ytop <- y + hicPlot$resolution

        if (half == "top") {
            hic_triangle <- polygonGrob(
                x = c(xleft, xleft, xright),
                y = c(ybottom, ytop, ytop),
                gp = gpar(col = NA, fill = col),
                default.units = "native"
            )
        } else if (half == "bottom") {
            hic_triangle <- polygonGrob(
                x = c(xleft, xright, xright),
                y = c(ybottom, ybottom, ytop),
                gp = gpar(col = NA, fill = col),
                default.units = "native"
            )
        }

        assign("hic_grobs",
            addGrob(
                gTree = get("hic_grobs", envir = pgEnv),
                child = hic_triangle
            ),
            envir = pgEnv
        )
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    hicInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "hicInternal"
    )
    
    ## Justification
    hicInternal$just <- justConversion(just = hicInternal$just)

    if (is.null(hicInternal$quiet)) hicInternal$quiet <- FALSE

    if (is.null(hicInternal$data)) stop("argument \"data\" is missing, ",
                                        "with no default.", call. = FALSE)
    if (is.null(hicInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    hicPlot <- structure(list(
        chrom = hicInternal$chrom,
        chromstart = hicInternal$chromstart,
        chromend = hicInternal$chromend,
        altchrom = hicInternal$altchrom,
        altchromstart = hicInternal$altchromstart,
        altchromend = hicInternal$altchromend,
        assembly = hicInternal$assembly,
        resolution = hicInternal$resolution,
        x = hicInternal$x, y = hicInternal$y,
        width = hicInternal$width,
        height = hicInternal$height,
        just = hicInternal$just,
        color_palette = NULL,
        colorTrans = hicInternal$colorTrans,
        zrange = hicInternal$zrange,
        half = hicInternal$half,
        grobs = NULL
    ), class = "hicSquare")
    attr(x = hicPlot, which = "plotted") <- hicInternal$draw

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    hicPlot$assembly <- parseAssembly(assembly = hicPlot$assembly)

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_placement(object = hicPlot)
    errorcheck_plotHicSquare(
        hic = hicInternal$data,
        hicPlot = hicPlot,
        norm = hicInternal$norm
    )

    # =========================================================================
    # PARSE UNITS AND Y-COORD
    # =========================================================================

    hicPlot <- defaultUnits(
        object = hicPlot,
        default.units = hicInternal$default.units
    )
    
    if (!is.null(hicPlot$x) & !is.null(hicPlot$y)){
        
        ## Convert width and height to same unit to check that they're the same
        width <- convertWidth(hicPlot$width, unitTo = get("page_units", 
                                                        envir = pgEnv))
        height <- convertHeight(hicPlot$height, unitTo = get("page_units",
                                                            envir = pgEnv))
        
        if (!identical(width, height)){
            stop("Attempting to plot square Hi-C plot with different width and",
                " height. Use `plotHicTriangle` or `plotHicRectangle` ",
                " for non-square Hi-C plots.", call. = FALSE)
        }
        
    }
    
    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================


    if (is.null(hicPlot$chromstart) & is.null(hicPlot$chromend)) {
        chromData <- get_wholeChrom(
            chrom = hicPlot$chrom,
            assembly = hicPlot$assembly
        )
        hicPlot$chromstart <- chromData[[1]]
        hicPlot$chromend <- chromData[[2]]
    }

    if (!is.null(hicPlot$altchrom)) {
        if (is.null(hicPlot$altchromstart) & is.null(hicPlot$altchromend)) {
            altchromData <- get_wholeChrom(
                chrom = hicPlot$altchrom,
                assembly = hicPlot$assembly
            )
            hicPlot$altchromstart <- altchromData[[1]]
            hicPlot$altchromend <- altchromData[[2]]
        }
    }

    # =========================================================================
    # ADJUST RESOLUTION
    # =========================================================================

    if (hicPlot$resolution == "auto") {
        hicPlot <- adjust_resolution(
            hic = hicInternal$data,
            hicPlot = hicPlot
        )
    } else {
        hicPlot <- hic_limit(
            hic = hicInternal$data,
            hicPlot = hicPlot
        )
    }

    # =========================================================================
    # READ IN DATA
    # =========================================================================

    hic <- read_data(
        hic = hicInternal$data,
        hicPlot = hicPlot,
        norm = hicInternal$norm,
        assembly = hicPlot$assembly,
        type = hicInternal$matrix,
        quiet = hicInternal$quiet
    )

    # =========================================================================
    # SUBSET DATA
    # =========================================================================

    hic <- subset_data(hic = hic, hicPlot = hicPlot)

    # =========================================================================
    # MAKE SYMMETRIC
    # =========================================================================

    hicFlip <- hic[, c("y", "x", "counts")]
    colnames(hicFlip) <- c("x", "y", "counts")
    hic <- unique(rbind(hic, hicFlip))
    colnames(hic) <- c("x", "y", "counts")

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

        ## Log color scale
        if (grepl("log", hicInternal$colorTrans) == TRUE) {
            logBase <- utils::type.convert(gsub("log", "", 
                                                hicInternal$colorTrans), 
                                as.is = TRUE)
            if (is.na(logBase)) {
                logBase <- exp(1)
            }

            ## Won't scale to log if negative values
            if (any(hic$counts < 0)) {
                stop("Negative values in Hi-C data. Cannot scale colors on a ",
                "log scale. Please set `colorTrans = 'linear'`.", call. = FALSE)
            }

            hic$counts <- log(hic$counts, base = logBase)
            hic$color <- mapColors(vector = hic$counts,
                palette = hicInternal$palette,
                range = log(hicPlot$zrange, logBase)
            )
        } else {
            hic$color <- mapColors(vector = hic$counts,
                palette = hicInternal$palette,
                range = hicPlot$zrange
            )
        }


        hicPlot$color_palette <- hicInternal$palette
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport xscale and yscale
    scale <- vp_scale(hicPlot = hicPlot)

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "hicSquare",
        length(grep(
            pattern = "hicSquare",
            x = currentViewports
        )) + 1
    )

    ## If placing information is provided but plot == TRUE, set up it's own
    ## viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(hicPlot$x) | is.null(hicPlot$y)) {
        vp <- viewport(
            height = unit(1, "snpc"), width = unit(1, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            clip = "on",
            xscale = scale[[1]], yscale = scale[[2]],
            just = "center",
            name = vp_name
        )

        if (hicInternal$draw == TRUE) {
            vp$name <- "hicSquare1"
            grid.newpage()
        }
    } else {
        addViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = hicPlot)

        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = scale[[1]], yscale = scale[[2]],
            just = hicInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS
    # =========================================================================

    assign("hic_grobs", gTree(vp = vp), envir = pgEnv)

    # =========================================================================
    # MAKE GROBS
    # =========================================================================

    if (!is.null(hicPlot$chromstart) & !is.null(hicPlot$chromend)) {
        if (nrow(hic) > 0) {

            ## Determine which grobs will be squares [[1]] and which
            ## will be triangles [[2]]
            shapes <- hic_shapes(
                hic = hic, hicPlot = hicPlot,
                half = hicInternal$half
            )

            if (!is.null(shapes[[1]])) {

                ## Make square grobs and add to grob gTree
                hic_squares <- rectGrob(
                    x = shapes[[1]]$x,
                    y = shapes[[1]]$y,
                    just = c("left", "bottom"),
                    width = hicPlot$resolution,
                    height = hicPlot$resolution,
                    gp = gpar(col = NA, fill = shapes[[1]]$color),
                    default.units = "native"
                )

                assign("hic_grobs",
                    addGrob(
                        gTree = get("hic_grobs", envir = pgEnv),
                        child = hic_squares
                    ),
                    envir = pgEnv
                )
            }

            ## Make triangle grobs and add to grob gTree
            if (!is.null(shapes[[2]])) {
                invisible(apply(shapes[[2]], 1, hic_diagonal,
                    hicPlot = hicPlot,
                    half = hicInternal$half
                ))
            }
        } else {
            
            if (!is.na(hicPlot$resolution)){
                warning("No data found in region. Suggestions: check that ",
                        "chromosome names match genome assembly; ",
                        "check region.", call. = FALSE)
            }
            
        }
    }


    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (hicInternal$draw == TRUE) {
        grid.draw(get("hic_grobs", envir = pgEnv))
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    hicPlot$grobs <- get("hic_grobs", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================
    if (is.null(hicPlot$altchrom)) {
        hicPlot$altchrom <- hicPlot$chrom
        hicPlot$altchromstart <- hicPlot$chromstart
        hicPlot$altchromend <- hicPlot$chromend
    }

    message("hicSquare[", vp$name, "]")
    invisible(hicPlot)
}
