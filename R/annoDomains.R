#' Annotate domains in a Hi-C plot
#'
#' @usage annoDomains(
#'     plot,
#'     data,
#'     half = "inherit",
#'     linecolor = "black",
#'     params = NULL,
#'     ...
#' )
#'
#' @param plot Hi-C plot object from \code{plotHicSquare} or
#' \code{plotHicTriangle} on which to annotate pixels.
#' @param data A string specifying the BED file path, a dataframe in BED
#' format, or a \link[GenomicRanges]{GRanges} object specifying
#' domain ranges.
#' @param half Character value specifying which half of hic plots
#' to annotate. Triangle Hi-C plots will always default to the entirety of
#' the triangular plot. Default value is \code{half = "inherit"}. Options are:
#' \itemize{
#' \item{\code{"inherit"}: }{Domains will be annotated on the \code{half}
#' inherited by the input Hi-C plot.}
#' \item{\code{"both"}: }{Domains will be annotated on both halves of the
#' diagonal of a square Hi-C plot.}
#' \item{\code{"top"}: }{Domains will be annotated on the upper diagonal
#' half of a square Hi-C plot.}
#' \item{\code{"bottom"}: }{Domains will be annotated ont the bottom diagonal
#' half of a square Hi-C plot.}
#' }
#' @param linecolor A character value specifying the color of the domain
#' annotations. Default value is \code{linecolor = "black"}.
#' @param params An optional \link[plotgardener]{params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{domain} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Define a GRanges object with TAD ranges
#' library(GenomicRanges)
#' library(IRanges)
#' domains <- GRanges("chr21",
#'     ranges = IRanges(
#'         start = c(28210000, 29085000, 29430000, 29700000),
#'         end = c(29085000, 29430000, 29700000, 30125000)
#'     )
#' )
#'
#' ## Load Hi-C data
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
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
#' ## Annotate domains on bottom half 0f Hi-C plot
#' annoDomains(
#'     plot = hicPlot, data = domains,
#'     half = "bottom", linecolor = "red"
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
annoDomains <- function(plot, data, half = "inherit",
                        linecolor = "black", params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that checks errors for annoDomains
    errorcheck_annoDomains <- function(hic, half) {

        ###### hic #####

        ## check type of input for hic
        if (!class(hic) %in% c(
            "hicSquare", "hicTriangle",
            "hicRectangle"
        )) {
            stop("Input plot must be a plot of class \'hicSquare\', ",
                "\'hicTriangle\', or \'hicRectangle\'.",
                call. = FALSE
            )
        }

        ###### half #####

        ## half needs to be a valid option
        if (!half %in% c("inherit", "both", "top", "bottom")) {
            stop("Invalid \'half\'.  Options are \'inherit\', ",
                "\'both\', \'top\', or \'bottom\'.",
                call. = FALSE
            )
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
                stop("Cannot annotate domains for an intrachromosomal ",
                    "square Hi-C plot.",
                    call. = FALSE
                )
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
    }

    ## Define a function that creates annotation grobs
    domainAnnotation <- function(df, object, half) {
        if (half == "bottom") {
            x0_1 <- utils::type.convert(df["start"], as.is = TRUE)
            x1_1 <- utils::type.convert(df["end"], as.is = TRUE)
            y_1 <- utils::type.convert(df["start"], as.is = TRUE)

            x_2 <- utils::type.convert(df["end"], as.is = TRUE)
            y0_2 <- utils::type.convert(df["start"], as.is = TRUE)
            y1_2 <- utils::type.convert(df["end"], as.is = TRUE)

            seg1 <- segmentsGrob(
                x0 = x0_1, x1 = x1_1,
                y0 = y_1, y1 = y_1,
                default.units = "native",
                gp = object$gp
            )

            seg2 <- segmentsGrob(
                x0 = x_2, x1 = x_2,
                y0 = y0_2, y1 = y1_2,
                default.units = "native",
                gp = object$gp
            )

            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = pgEnv),
                    child = seg1
                ),
                envir = pgEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = pgEnv),
                    child = seg2
                ),
                envir = pgEnv
            )
        } else if (half == "top") {
            x_1 <- utils::type.convert(df["start"], as.is = TRUE)
            y0_1 <- utils::type.convert(df["start"], as.is = TRUE)
            y1_1 <- utils::type.convert(df["end"], as.is = TRUE)

            x0_2 <- utils::type.convert(df["start"], as.is = TRUE)
            x1_2 <- utils::type.convert(df["end"], as.is = TRUE)
            y_2 <- utils::type.convert(df["end"], as.is = TRUE)

            seg1 <- segmentsGrob(
                x0 = x_1, x1 = x_1,
                y0 = y0_1, y1 = y1_1,
                default.units = "native",
                gp = object$gp
            )

            seg2 <- segmentsGrob(
                x0 = x0_2, x1 = x1_2,
                y0 = y_2, y1 = y_2,
                default.units = "native",
                gp = object$gp
            )

            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = pgEnv),
                    child = seg1
                ),
                envir = pgEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = pgEnv),
                    child = seg2
                ),
                envir = pgEnv
            )
        } else if (half == "both") {
            x0_1 <- utils::type.convert(df["start"], as.is = TRUE)
            x1_1 <- utils::type.convert(df["end"], as.is = TRUE)
            y_1 <- utils::type.convert(df["start"], as.is = TRUE)
            x_2 <- utils::type.convert(df["end"], as.is = TRUE)
            y0_2 <- utils::type.convert(df["start"], as.is = TRUE)
            y1_2 <- utils::type.convert(df["end"], as.is = TRUE)


            x_1 <- utils::type.convert(df["start"], as.is = TRUE)
            y0_1 <- utils::type.convert(df["start"], as.is = TRUE)
            y1_1 <- utils::type.convert(df["end"], as.is = TRUE)
            x0_2 <- utils::type.convert(df["start"], as.is = TRUE)
            x1_2 <- utils::type.convert(df["end"], as.is = TRUE)
            y_2 <- utils::type.convert(df["end"], as.is = TRUE)


            seg1 <- segmentsGrob(
                x0 = x0_1, x1 = x1_1,
                y0 = y_1, y1 = y_1,
                default.units = "native",
                gp = object$gp
            )

            seg2 <- segmentsGrob(
                x0 = x_2, x1 = x_2,
                y0 = y0_2, y1 = y1_2,
                default.units = "native",
                gp = object$gp
            )

            seg3 <- segmentsGrob(
                x0 = x_1, x1 = x_1,
                y0 = y0_1, y1 = y1_1,
                default.units = "native",
                gp = object$gp
            )

            seg4 <- segmentsGrob(
                x0 = x0_2, x1 = x1_2,
                y0 = y_2, y1 = y_2,
                default.units = "native",
                gp = object$gp
            )

            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = pgEnv),
                    child = seg1
                ),
                envir = pgEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = pgEnv),
                    child = seg2
                ),
                envir = pgEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = pgEnv),
                    child = seg3
                ),
                envir = pgEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = pgEnv),
                    child = seg4
                ),
                envir = pgEnv
            )
        }
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    domainsInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(
            match.call()[[1]]
        )),
        declaredArgs = lapply(
            match.call()[-1], eval.parent,
            n = 2
        ),
        class = "domainsInternal"
    )
    domainsInternal$gp <- setGP(
        gpList = gpar(),
        params = domainsInternal, ...
    )

    domainsInternal$gp$col <- domainsInternal$linecolor
    domainsInternal$gp$lineend <- "square"

    # =========================================================================
    # INITIALIZE OBJECT: GET REGION/DIMENSIONS FROM HIC PLOT INPUT
    # =========================================================================

    domains <- structure(list(
        hicClass = class(domainsInternal$plot),
        chrom = domainsInternal$plot$chrom,
        chromstart = domainsInternal$plot$chromstart,
        chromend = domainsInternal$plot$chromend,
        assembly = domainsInternal$plot$assembly,
        x = domainsInternal$plot$x,
        y = domainsInternal$plot$y,
        width = domainsInternal$plot$width,
        height = domainsInternal$plot$height,
        just = domainsInternal$plot$just, grobs = NULL,
        outsideVP = NULL
    ),
    class = "domain"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_page(error = "Cannot annotate Hi-C domains without a
                `plotgardener` page.")
    if (is.null(domainsInternal$plot)) {
        stop("argument \"plot\" is missing, ",
            "with no default.",
            call. = FALSE
        )
    }
    if (is.null(domainsInternal$data)) {
        stop("argument \"data\" is missing, ",
            "with no default.",
            call. = FALSE
        )
    }

    errorcheck_annoDomains(
        hic = domainsInternal$plot,
        half = domainsInternal$half
    )

    # =========================================================================
    # PARSE INHERITED HALF
    # =========================================================================

    half <- domainsInternal$half
    if (half == "inherit") {
        half <- inherit_half(hic = domainsInternal$plot)
    }

    if (is(domainsInternal$plot, "hicTriangle")  |
        is(domainsInternal$plot, "hicRectangle")) {
        half <- "top"
    }

    # =========================================================================
    # READ IN FILE, DATAFRAME OR GRANGES
    # =========================================================================

    bed <- read_rangeData(
        data = domainsInternal$data,
        assembly = domains$assembly
    )
    
    ## chrom format and data chrom format
    chromDataAgreement(data = bed, chrom = domains$chrom,
                    type = "ranges")

    # =========================================================================
    # SUBSET FOR DOMAINS IN REGION
    # =========================================================================

    bed_subset <- bed[which(bed[, "chrom"] == domains$chrom 
                            & bed[, "start"] <=
        domains$chromend & bed[, "end"] >=
        domains$chromstart), ]

    # =========================================================================
    # VIEWPORTS
    # =========================================================================
    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "domain",
        length(grep(
            pattern = "domain",
            x = currentViewports
        )) + 1
    )
    ## Make viewport based on hic input viewport
    if (is(domainsInternal$plot, "hicSquare")) {
        vp <- viewport(
            height = domainsInternal$plot$grobs$vp$height,
            width = domainsInternal$plot$grobs$vp$width,
            x = domainsInternal$plot$grobs$vp$x,
            y = domainsInternal$plot$grobs$vp$y,
            clip = "on",
            xscale = domainsInternal$plot$grobs$vp$xscale,
            yscale = domainsInternal$plot$grobs$vp$yscale,
            just = domainsInternal$plot$grobs$vp$justification,
            name = paste0(vp_name, "_inside")
        )

        vpClip <- NULL
    } else if (is(domainsInternal$plot, "hicTriangle")) {
        width <- convertUnit(domainsInternal$plot$outsideVP$width,
            unitTo = get("page_units", pgEnv), valueOnly = TRUE
        )

        vp <- viewport(
            height = unit(width / sqrt(2), get("page_units", pgEnv)),
            width = unit(width / sqrt(2), get("page_units", pgEnv)),
            x = unit(0, "npc"), y = unit(0, "npc"),
            xscale = domainsInternal$plot$grobs$vp$xscale,
            yscale = domainsInternal$plot$grobs$vp$yscale,
            just = c("left", "bottom"),
            name = paste0(vp_name, "_inside"),
            angle = -45
        )

        vpClip <- viewport(
            height = domainsInternal$plot$outsideVP$height,
            width = domainsInternal$plot$outsideVP$width,
            x = domainsInternal$plot$outsideVP$x,
            y = domainsInternal$plot$outsideVP$y,
            just = domainsInternal$plot$outsideVP$justification,
            clip = "on",
            name = paste0(vp_name, "_outside")
        )
    } else if (is(domainsInternal$plot, "hicRectangle")) {
        side <- convertUnit(domainsInternal$plot$grobs$vp$width,
            unitTo = get("page_units", pgEnv)
        )

        vp <- viewport(
            height = side, width = side,
            x = unit(0, "npc"), y = unit(0, "npc"),
            xscale = domainsInternal$plot$grobs$vp$xscale,
            yscale = domainsInternal$plot$grobs$vp$yscale,
            just = c("left", "bottom"),
            name = paste0(vp_name, "_inside"),
            angle = -45
        )

        vpClip <- viewport(
            height = domainsInternal$plot$outsideVP$height,
            width = domainsInternal$plot$outsideVP$width,
            x = domainsInternal$plot$outsideVP$x,
            y = domainsInternal$plot$outsideVP$y,
            just = domainsInternal$plot$outsideVP$justification,
            clip = "on",
            name = paste0(vp_name, "_outside")
        )
    }

    # =========================================================================
    # INITIALIZE GTREE
    # =========================================================================

    domains$outsideVP <- vpClip
    assign("domain_grobs", gTree(vp = vp), envir = pgEnv)

    # =========================================================================
    # GROBS
    # =========================================================================

    if (nrow(bed_subset) > 0) {
        invisible(apply(bed_subset, 1, domainAnnotation,
            object = domainsInternal,
            half = half
        ))
    } else {
        warning("No domains found in region.", call. = FALSE)
    }

    # =========================================================================
    # ADD GROBS TO GLIST AND OBJECT
    # =========================================================================

    domains$grobs <- get("domain_grobs", envir = pgEnv)

    if (!is.null(vpClip)) {
        pushViewport(vpClip)
        grid.draw(domains$grobs)
        upViewport()
    } else {
        grid.draw(domains$grobs)
    }


    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("domain[", vp_name, "]")
    invisible(domains)
}
