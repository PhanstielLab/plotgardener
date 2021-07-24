#' Annotate domains in a Hi-C plot
#'
#' @usage bb_annoDomains(
#'     plot,
#'     data,
#'     half = "inherit",
#'     linecolor = "black",
#'     params = NULL,
#'     ...
#' )
#'
#' @param plot Hi-C plot object from \code{bb_plotHicSquare} or
#' \code{bb_plotHicTriangle} on which to annotate pixels.
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
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_domain} object containing relevant
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
#' library(BentoBoxData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 4.5, height = 4, default.units = "inches")
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- bb_plotHicSquare(
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
#' bb_annoDomains(
#'     plot = hicPlot, data = domains,
#'     half = "bottom", linecolor = "red"
#' )
#'
#' ## Annotate heatmap legend
#' bb_annoHeatmapLegend(
#'     plot = hicPlot,
#'     x = 3.6, y = 0.5, width = 0.12, height = 1.2,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bb_annoGenomeLabel(
#'     plot = hicPlot, x = 0.5, y = 3.53, scale = "Mb",
#'     just = c("left", "top")
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @export
bb_annoDomains <- function(plot, data, half = "inherit",
                        linecolor = "black", params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that checks errors for bb_annoDomains
    errorcheck_bb_annoDomains <- function(hic, half) {

        ###### hic #####

        ## check type of input for hic
        if (!class(hic) %in% c(
            "bb_hicSquare", "bb_hicTriangle",
            "bb_hicRectangle"
        )) {
            stop("Input plot must be a plot of class \'bb_hicSquare\', ",
                "\'bb_hicTriangle\', or \'bb_hicRectangle\'.",
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
                stop("Cannot annotate domains for an intrachromosomal ",
                    "square Hi-C plot.",
                    call. = FALSE
                )
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
    }

    ## Define a function that creates annotatino grobs
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
                    get("domain_grobs", envir = bbEnv),
                    child = seg1
                ),
                envir = bbEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = bbEnv),
                    child = seg2
                ),
                envir = bbEnv
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
                    get("domain_grobs", envir = bbEnv),
                    child = seg1
                ),
                envir = bbEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = bbEnv),
                    child = seg2
                ),
                envir = bbEnv
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
                    get("domain_grobs", envir = bbEnv),
                    child = seg1
                ),
                envir = bbEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = bbEnv),
                    child = seg2
                ),
                envir = bbEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = bbEnv),
                    child = seg3
                ),
                envir = bbEnv
            )
            assign("domain_grobs",
                addGrob(
                    get("domain_grobs", envir = bbEnv),
                    child = seg4
                ),
                envir = bbEnv
            )
        }
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_domainsInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(
            match.call()[[1]]
        )),
        declaredArgs = lapply(
            match.call()[-1], eval.parent,
            n = 2
        ),
        class = "bb_domainsInternal"
    )
    bb_domainsInternal$gp <- setGP(
        gpList = gpar(),
        params = bb_domainsInternal, ...
    )

    bb_domainsInternal$gp$col <- bb_domainsInternal$linecolor
    bb_domainsInternal$gp$lineend <- "square"

    # =========================================================================
    # INITIALIZE OBJECT: GET REGION/DIMENSIONS FROM HIC PLOT INPUT
    # =========================================================================

    bb_domains <- structure(list(
        hicClass = class(bb_domainsInternal$plot),
        chrom = bb_domainsInternal$plot$chrom,
        chromstart = bb_domainsInternal$plot$chromstart,
        chromend = bb_domainsInternal$plot$chromend,
        assembly = bb_domainsInternal$plot$assembly,
        x = bb_domainsInternal$plot$x,
        y = bb_domainsInternal$plot$y,
        width = bb_domainsInternal$plot$width,
        height = bb_domainsInternal$plot$height,
        just = bb_domainsInternal$plot$just, grobs = NULL,
        outsideVP = NULL
    ),
    class = "bb_domain"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_bbpage(error = "Cannot annotate Hi-C domains without a
                BentoBox page.")
    if (is.null(bb_domainsInternal$plot)) {
        stop("argument \"plot\" is missing, ",
            "with no default.",
            call. = FALSE
        )
    }
    if (is.null(bb_domainsInternal$data)) {
        stop("argument \"data\" is missing, ",
            "with no default.",
            call. = FALSE
        )
    }

    errorcheck_bb_annoDomains(
        hic = bb_domainsInternal$plot,
        half = bb_domainsInternal$half
    )

    # =========================================================================
    # PARSE INHERITED HALF
    # =========================================================================

    half <- bb_domainsInternal$half
    if (half == "inherit") {
        half <- inherit_half(hic = bb_domainsInternal$plot)
    }

    if (is(bb_domainsInternal$plot, "bb_hicTriangle")  |
        is(bb_domainsInternal$plot, "bb_hicRectangle")) {
        half <- "top"
    }

    # =========================================================================
    # READ IN FILE, DATAFRAME OR GRANGES
    # =========================================================================

    bed <- read_rangeData(
        data = bb_domainsInternal$data,
        assembly = bb_domains$assembly
    )
    
    ## chrom format and data chrom format
    chromDataAgreement(data = bed, chrom = bb_domains$chrom,
                    type = "ranges")

    # =========================================================================
    # SUBSET FOR DOMAINS IN REGION
    # =========================================================================

    bed_subset <- bed[which(bed[, "chrom"] == bb_domains$chrom 
                            & bed[, "start"] <=
        bb_domains$chromend & bed[, "end"] >=
        bb_domains$chromstart), ]

    # =========================================================================
    # VIEWPORTS
    # =========================================================================
    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_domain",
        length(grep(
            pattern = "bb_domain",
            x = currentViewports
        )) + 1
    )
    ## Make viewport based on hic input viewport
    if (is(bb_domainsInternal$plot, "bb_hicSquare")) {
        vp <- viewport(
            height = bb_domainsInternal$plot$grobs$vp$height,
            width = bb_domainsInternal$plot$grobs$vp$width,
            x = bb_domainsInternal$plot$grobs$vp$x,
            y = bb_domainsInternal$plot$grobs$vp$y,
            clip = "on",
            xscale = bb_domainsInternal$plot$grobs$vp$xscale,
            yscale = bb_domainsInternal$plot$grobs$vp$yscale,
            just = bb_domainsInternal$plot$grobs$vp$justification,
            name = paste0(vp_name, "_inside")
        )

        vpClip <- NULL
    } else if (is(bb_domainsInternal$plot, "bb_hicTriangle")) {
        width <- convertUnit(bb_domainsInternal$plot$outsideVP$width,
            unitTo = get("page_units", bbEnv), valueOnly = TRUE
        )

        vp <- viewport(
            height = unit(width / sqrt(2), get("page_units", bbEnv)),
            width = unit(width / sqrt(2), get("page_units", bbEnv)),
            x = unit(0, "npc"), y = unit(0, "npc"),
            xscale = bb_domainsInternal$plot$grobs$vp$xscale,
            yscale = bb_domainsInternal$plot$grobs$vp$yscale,
            just = c("left", "bottom"),
            name = paste0(vp_name, "_inside"),
            angle = -45
        )

        vpClip <- viewport(
            height = bb_domainsInternal$plot$outsideVP$height,
            width = bb_domainsInternal$plot$outsideVP$width,
            x = bb_domainsInternal$plot$outsideVP$x,
            y = bb_domainsInternal$plot$outsideVP$y,
            just = bb_domainsInternal$plot$outsideVP$justification,
            clip = "on",
            name = paste0(vp_name, "_outside")
        )
    } else if (is(bb_domainsInternal$plot, "bb_hicRectangle")) {
        side <- convertUnit(bb_domainsInternal$plot$grobs$vp$width,
            unitTo = get("page_units", bbEnv)
        )

        vp <- viewport(
            height = side, width = side,
            x = unit(0, "npc"), y = unit(0, "npc"),
            xscale = bb_domainsInternal$plot$grobs$vp$xscale,
            yscale = bb_domainsInternal$plot$grobs$vp$yscale,
            just = c("left", "bottom"),
            name = paste0(vp_name, "_inside"),
            angle = -45
        )

        vpClip <- viewport(
            height = bb_domainsInternal$plot$outsideVP$height,
            width = bb_domainsInternal$plot$outsideVP$width,
            x = bb_domainsInternal$plot$outsideVP$x,
            y = bb_domainsInternal$plot$outsideVP$y,
            just = bb_domainsInternal$plot$outsideVP$justification,
            clip = "on",
            name = paste0(vp_name, "_outside")
        )
    }

    # =========================================================================
    # INITIALIZE GTREE
    # =========================================================================

    bb_domains$outsideVP <- vpClip
    assign("domain_grobs", gTree(vp = vp), envir = bbEnv)

    # =========================================================================
    # GROBS
    # =========================================================================

    if (nrow(bed_subset) > 0) {
        invisible(apply(bed_subset, 1, domainAnnotation,
            object = bb_domainsInternal,
            half = half
        ))
    } else {
        warning("No domains found in region.", call. = FALSE)
    }

    # =========================================================================
    # ADD GROBS TO GLIST AND OBJECT
    # =========================================================================

    bb_domains$grobs <- get("domain_grobs", envir = bbEnv)

    if (!is.null(vpClip)) {
        pushViewport(vpClip)
        grid.draw(bb_domains$grobs)
        upViewport()
    } else {
        grid.draw(bb_domains$grobs)
    }


    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_domain[", vp_name, "]")
    invisible(bb_domains)
}
