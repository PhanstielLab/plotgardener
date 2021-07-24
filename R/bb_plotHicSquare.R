#' Plot a Hi-C interaction matrix in a square format
#' 
#' @usage bb_plotHicSquare(
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
#' \link[BentoBox]{bb_assembly} object.
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
#' plotted BentoBox plot according to the units of the BentoBox page.
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
#' @param params An optional \link[BentoBox]{bb_params} object containing
#' relevant function parameters.
#' @param quiet A logical indicating whether or not to print messages.
#'
#' @return Returns a \code{bb_hicSquare} object containing relevant
#' genomic region, Hi-C data, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data
#' library(BentoBoxData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create a page
#' bb_pageCreate(width = 3, height = 3, default.units = "inches")
#'
#' ## Plot and place Hi-C plot
#' hicPlot <- bb_plotHicSquare(
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
#' bb_annoHeatmapLegend(
#'     plot = hicPlot, x = 2.6, y = 0.5,
#'     width = 0.12, height = 1.2,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate x-axis and y-axis genome labels
#' bb_annoGenomeLabel(
#'     plot = hicPlot, scale = "Mb", axis = "x",
#'     x = 0.5, y = 2.53, just = c("left", "top")
#' )
#' bb_annoGenomeLabel(
#'     plot = hicPlot, scale = "Mb", axis = "y",
#'     x = 0.47, y = 0.5, just = c("right", "top")
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @details
#' A square Hi-C plot can be placed on a BentoBox coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' bb_plotHicSquare(data, chrom,
#'                 chromstart = NULL, chromend = NULL,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can be used to quickly plot an unannotated
#' square Hi-C plot by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotHicSquare(data, chrom,
#'                 chromstart = NULL, chromend = NULL)
#' }
#'
#' @seealso \link[BentoBox]{bb_readHic}
#'
#' @export
bb_plotHicSquare <- function(data, resolution = "auto", zrange = NULL,
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

    ## Define a function that catches errors for bb_plothic
    errorcheck_bb_plothic <- function(hic, hic_plot, norm) {
        
        ###### hic/norm 
        bb_hicErrors(hic = hic, norm = norm)

        ###### chrom/chromstart/chromend/
        ###### altchrom/altchromstart/altchromend

        ## Even though straw technically works without "chr" for hg19,
        ## will not accept for consistency purposes
        if (hic_plot$assembly$Genome == "hg19") {
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
        
        if (!is.null(hic_plot$altchrom)) {

            ## Can't specify altchrom without a chrom
            if (is.null(hic_plot$chrom)) {
                stop("Specified \'altchrom\', but did not give \'chrom\'.",
                    call. = FALSE
                )
            }
            ## Even though straw technically works without "chr" for hg19,
            ## will not accept for consistency purposes
            if (hic_plot$assembly$Genome == "hg19") {
                if (grepl("chr", hic_plot$altchrom) == FALSE) {
                    stop("'", hic_plot$altchrom, "'",
                        "is an invalid input for an hg19 chromsome. ",
                        "Please specify chromosome as",
                        "'chr", hic_plot$altchrom, "'.",
                        call. = FALSE
                    )
                }
            }

            bb_regionErrors(chromstart = hic_plot$altchromstart,
                        chromend = hic_plot$altchromend)


            if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend) &
                !is.null(hic_plot$altchromstart) &
                !is.null(hic_plot$altchromend)) {

                ## Check to see if region is square
                if ((hic_plot$chromend - hic_plot$chromstart) !=
                    (hic_plot$altchromend - hic_plot$altchromstart)) {
                    warning("Trying to plot non-square region.", call. = FALSE)
                }
            }
        }

        ###### zrange #####
        bb_rangeErrors(range = hic_plot$zrange)

        ###### half #####

        if (is.null(hic_plot$altchrom)) {
            if (!(hic_plot$half %in% c("both", "top", "bottom"))) {
                stop("Invalid \'half\'.  Options are \'both\', ",
                    "\top\', or \'bottom\'.",
                    call. = FALSE
                )
            }
        } else {
            if (!hic_plot$half %in% c("top", "bottom")) {
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
    subset_data <- function(hic, hic_plot) {
        if (nrow(hic) > 0) {
            if (is.null(hic_plot$altchrom)) {
                hic <- hic[which(hic[, "x"] >= hic_plot$chromstart -
                    hic_plot$resolution &
                    hic[, "x"] <= hic_plot$chromend + hic_plot$resolution &
                    hic[, "y"] >= hic_plot$chromstart - hic_plot$resolution &
                    hic[, "y"] <= hic_plot$chromend + hic_plot$resolution), ]
            } else {
                if (hic_plot$chrom != hic_plot$altchrom) {
                    hic <- hic[which(hic[, "x"] >= hic_plot$chromstart -
                        hic_plot$resolution &
                        hic[, "x"] <= hic_plot$chromend + hic_plot$resolution &
                        hic[, "y"] >= hic_plot$altchromstart -
                            hic_plot$resolution &
                        hic[, "y"] <= hic_plot$altchromend +
                            hic_plot$resolution), ]
                }
            }
        }

        return(hic)
    }

    ## Define a function that sets viewport xscale and yscale
    vp_scale <- function(hic_plot) {
        if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)) {
            if (is.null(hic_plot$altchrom)) {
                xscale <- c(hic_plot$chromstart, hic_plot$chromend)
                yscale <- xscale
            } else {
                if (hic_plot$half == "bottom") {
                    xscale <- c(hic_plot$chromstart, hic_plot$chromend)
                    yscale <- c(hic_plot$altchromstart, hic_plot$altchromend)
                }

                if (hic_plot$half == "top") {
                    xscale <- c(hic_plot$altchromstart, hic_plot$altchromend)
                    yscale <- c(hic_plot$chromstart, hic_plot$chromend)
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
    hic_shapes <- function(hic, hic_plot, half) {
        if (!is.null(hic_plot$altchrom)) {

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
    hic_diagonal <- function(hic, hic_plot, half) {
        col <- hic["color"]
        x <- utils::type.convert(hic["x"], as.is = TRUE)
        y <- utils::type.convert(hic["y"], as.is = TRUE)

        xleft <- x
        xright <- x + hic_plot$resolution
        ybottom <- y
        ytop <- y + hic_plot$resolution

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
                gTree = get("hic_grobs", envir = bbEnv),
                child = hic_triangle
            ),
            envir = bbEnv
        )
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_hicInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_hicInternal"
    )
    
    ## Justification
    bb_hicInternal$just <- bb_justConversion(just = bb_hicInternal$just)

    if (is.null(bb_hicInternal$quiet)) bb_hicInternal$quiet <- FALSE

    if (is.null(bb_hicInternal$data)) stop("argument \"data\" is missing, ",
                                        "with no default.", call. = FALSE)
    if (is.null(bb_hicInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    hic_plot <- structure(list(
        chrom = bb_hicInternal$chrom,
        chromstart = bb_hicInternal$chromstart,
        chromend = bb_hicInternal$chromend,
        altchrom = bb_hicInternal$altchrom,
        altchromstart = bb_hicInternal$altchromstart,
        altchromend = bb_hicInternal$altchromend,
        assembly = bb_hicInternal$assembly,
        resolution = bb_hicInternal$resolution,
        x = bb_hicInternal$x, y = bb_hicInternal$y,
        width = bb_hicInternal$width,
        height = bb_hicInternal$height,
        just = bb_hicInternal$just,
        color_palette = NULL,
        colorTrans = bb_hicInternal$colorTrans,
        zrange = bb_hicInternal$zrange,
        half = bb_hicInternal$half,
        grobs = NULL
    ), class = "bb_hicSquare")
    attr(x = hic_plot, which = "plotted") <- bb_hicInternal$draw

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    hic_plot$assembly <- parse_bbAssembly(assembly = hic_plot$assembly)

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_placement(object = hic_plot)
    errorcheck_bb_plothic(
        hic = bb_hicInternal$data,
        hic_plot = hic_plot,
        norm = bb_hicInternal$norm
    )

    # =========================================================================
    # PARSE UNITS AND Y-COORD
    # =========================================================================

    hic_plot <- defaultUnits(
        object = hic_plot,
        default.units = bb_hicInternal$default.units
    )
    
    if (!is.null(hic_plot$x) & !is.null(hic_plot$y)){
        
        ## Convert width and height to same unit to check that they're the same
        width <- convertWidth(hic_plot$width, unitTo = get("page_units", 
                                                        envir = bbEnv))
        height <- convertHeight(hic_plot$height, unitTo = get("page_units",
                                                            envir = bbEnv))
        
        if (!identical(width, height)){
            stop("Attempting to plot square Hi-C plot with different width and",
                " height. Use `bb_plotHicTriangle` or `bb_plotHicRectangle` ",
                " for non-square Hi-C plots.", call. = FALSE)
        }
        
    }
    
    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================


    if (is.null(hic_plot$chromstart) & is.null(hic_plot$chromend)) {
        chromData <- get_wholeChrom(
            chrom = hic_plot$chrom,
            assembly = hic_plot$assembly
        )
        hic_plot$chromstart <- chromData[[1]]
        hic_plot$chromend <- chromData[[2]]
    }

    if (!is.null(hic_plot$altchrom)) {
        if (is.null(hic_plot$altchromstart) & is.null(hic_plot$altchromend)) {
            altchromData <- get_wholeChrom(
                chrom = hic_plot$altchrom,
                assembly = hic_plot$assembly
            )
            hic_plot$altchromstart <- altchromData[[1]]
            hic_plot$altchromend <- altchromData[[2]]
        }
    }

    # =========================================================================
    # ADJUST RESOLUTION
    # =========================================================================

    if (hic_plot$resolution == "auto") {
        hic_plot <- adjust_resolution(
            hic = bb_hicInternal$data,
            hic_plot = hic_plot
        )
    } else {
        hic_plot <- hic_limit(
            hic = bb_hicInternal$data,
            hic_plot = hic_plot
        )
    }

    # =========================================================================
    # READ IN DATA
    # =========================================================================

    hic <- read_data(
        hic = bb_hicInternal$data,
        hic_plot = hic_plot,
        norm = bb_hicInternal$norm,
        assembly = hic_plot$assembly,
        type = bb_hicInternal$matrix,
        quiet = bb_hicInternal$quiet
    )

    # =========================================================================
    # SUBSET DATA
    # =========================================================================

    hic <- subset_data(hic = hic, hic_plot = hic_plot)

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

    hic_plot <- set_zrange(hic = hic, hic_plot = hic_plot)
    hic$counts[hic$counts <= hic_plot$zrange[1]] <- hic_plot$zrange[1]
    hic$counts[hic$counts >= hic_plot$zrange[2]] <- hic_plot$zrange[2]

    # =========================================================================
    # CONVERT NUMBERS TO COLORS
    # =========================================================================

    ## if we don't have an appropriate zrange (even after setting it
    ## based on a null zrange), can't scale to colors
    if (!is.null(hic_plot$zrange) & length(unique(hic_plot$zrange)) == 2) {

        ## Log color scale
        if (grepl("log", bb_hicInternal$colorTrans) == TRUE) {
            logBase <- utils::type.convert(gsub("log", "", 
                                                bb_hicInternal$colorTrans), 
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
            hic$color <- bb_mapColors(vector = hic$counts,
                palette = bb_hicInternal$palette,
                range = log(hic_plot$zrange, logBase)
            )
        } else {
            hic$color <- bb_mapColors(vector = hic$counts,
                palette = bb_hicInternal$palette,
                range = hic_plot$zrange
            )
        }


        hic_plot$color_palette <- bb_hicInternal$palette
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport xscale and yscale
    scale <- vp_scale(hic_plot = hic_plot)

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_hicSquare",
        length(grep(
            pattern = "bb_hicSquare",
            x = currentViewports
        )) + 1
    )

    ## If placing information is provided but plot == TRUE, set up it's own
    ## viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(hic_plot$x) | is.null(hic_plot$y)) {
        vp <- viewport(
            height = unit(1, "snpc"), width = unit(1, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            clip = "on",
            xscale = scale[[1]], yscale = scale[[2]],
            just = "center",
            name = vp_name
        )

        if (bb_hicInternal$draw == TRUE) {
            vp$name <- "bb_hicSquare1"
            grid.newpage()
        }
    } else {
        add_bbViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = hic_plot)

        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = scale[[1]], yscale = scale[[2]],
            just = bb_hicInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS
    # =========================================================================

    assign("hic_grobs", gTree(vp = vp), envir = bbEnv)

    # =========================================================================
    # MAKE GROBS
    # =========================================================================

    if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)) {
        if (nrow(hic) > 0) {

            ## Determine which grobs will be squares [[1]] and which
            ## will be triangles [[2]]
            shapes <- hic_shapes(
                hic = hic, hic_plot = hic_plot,
                half = bb_hicInternal$half
            )

            if (!is.null(shapes[[1]])) {

                ## Make square grobs and add to grob gTree
                hic_squares <- rectGrob(
                    x = shapes[[1]]$x,
                    y = shapes[[1]]$y,
                    just = c("left", "bottom"),
                    width = hic_plot$resolution,
                    height = hic_plot$resolution,
                    gp = gpar(col = NA, fill = shapes[[1]]$color),
                    default.units = "native"
                )

                assign("hic_grobs",
                    addGrob(
                        gTree = get("hic_grobs", envir = bbEnv),
                        child = hic_squares
                    ),
                    envir = bbEnv
                )
            }

            ## Make triangle grobs and add to grob gTree
            if (!is.null(shapes[[2]])) {
                invisible(apply(shapes[[2]], 1, hic_diagonal,
                    hic_plot = hic_plot,
                    half = bb_hicInternal$half
                ))
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
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (bb_hicInternal$draw == TRUE) {
        grid.draw(get("hic_grobs", envir = bbEnv))
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    hic_plot$grobs <- get("hic_grobs", envir = bbEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================
    if (is.null(hic_plot$altchrom)) {
        hic_plot$altchrom <- hic_plot$chrom
        hic_plot$altchromstart <- hic_plot$chromstart
        hic_plot$altchromend <- hic_plot$chromend
    }

    message("bb_hicSquare[", vp$name, "]")
    invisible(hic_plot)
}
