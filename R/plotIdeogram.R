#' Plot a chromosome ideogram with or without cytobands
#' 
#' @usage plotIdeogram(
#'     chrom,
#'     assembly = "hg38",
#'     data = NULL,
#'     orientation = "h",
#'     showBands = TRUE,
#'     fill = NULL,
#'     x = NULL,
#'     y = NULL,
#'     width = NULL,
#'     height = NULL,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     draw = TRUE,
#'     params = NULL
#' )
#'
#' @param chrom Chromosome to be plotted, as a string.
#' @param assembly Default genome assembly as a string or a
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param data Custom cytoband data, as a dataframe with the following
#' columns: "seqnames", "start", "end", "width", "strand", 
#' "name", "gieStain".
#' @param orientation Character value indicating the orientation
#' of the ideogram. Default value is \code{orientation = "h"}.
#' Options are:
#' \itemize{
#' \item{\code{"v"}: }{Vertical ideogram orientation.}
#' \item{\code{"h"}: }{Horizontal ideogram orientation.}
#' }
#' @param showBands Logical value indicating whether to draw
#' colored cytobands within ideogram.
#' Default value is \code{showBands = TRUE}.
#' @param fill (optional) A vector specifying alternate colors for cytoband
#' stains. To change specific gieStain values (i.e. gneg, gpos, etc.) to 
#' specific colors, this vector can be named. This vector must have the same
#' number of colors as there are gieStain values for each genome assembly.
#' @param x A numeric or unit object specifying ideogram x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying ideogram y-location.
#' The character value will
#' place the ideogram y relative to the bottom of the most recently
#' plotted plot according to the units of the plotgardener page.
#' @param width A numeric or unit object specifying ideogram width.
#' @param height A numeric or unit object specifying ideogram height.
#' @param just Justification of ideogram relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal justification
#' and the second value specifies vertical justification.
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
#'
#' @return Returns a \code{ideogram} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load Giemsa stain band information and genomic
#' ## annotation data for hg19 genome assembly
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(AnnotationHub)
#'
#' ## Create page
#' pageCreate(width = 4.5, height = 1, default.units = "inches")
#'
#' ## Plot and place ideogram
#' ideogramPlot <- plotIdeogram(
#'     chrom = "chr2", assembly = "hg19",
#'     x = 0.25, y = 0.25, width = 4, height = 0.3,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Plot text
#' plotText(
#'     label = "Chromosome 2", fontcolor = "dark grey",
#'     x = 4.25, y = 0.65, just = "right"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' An ideogram can be placed on a plotgardener coordinate page by
#' providing plot placement parameters:
#' \preformatted{
#' plotIdeogram(chrom,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated ideogram
#' by ignoring plot placement parameters:
#' \preformatted{
#' plotIdeogram(chrom)
#' }
#' If no data is provided, Giemsa stain band data will first try to 
#' fetch UCSC with AnnotationHub. The results are cached for faster access, 
#' but these cached items can be deleted. If no internet connection is 
#' available and AnnotationHub has not previously cached the data, custom 
#' Giemsa stain band data can be loaded with the `data` parameter.
#'
#' @seealso \link[AnnotationHub]{AnnotationHub}
#' 
#' @export
plotIdeogram <- function(chrom, assembly = "hg38", data = NULL,
                            orientation = "h", showBands = TRUE, fill = NULL,
                            x = NULL, y = NULL,
                            width = NULL, height = NULL,
                            just = c("left", "top"), default.units = "inches",
                            draw = TRUE, params = NULL) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that checks errors for plotIdeogram
    errorcheck_plotIdeogram <- function(orientation) {
        if (!orientation %in% c("v", "h")) {
            stop("Invalid /'orientation/' parameter. Options are 'v' or 'h'.",
                call. = FALSE
            )
        }
    }

    ## Define a function to get cytoBand data for a genome assembly
    cytoAssembly <- function(assembly) {
    
        ## Get string assembly name
        assemblyName <- assembly$Genome
        
        ## Check in defaults for data
        if(!any(cytoband_AH_assembly$assembly %in% assemblyName)){
            warning("UCSC cytoBand data not available for the given genome",
                    " assembly. Default data can only be found for the", 
                    " following assemblies:",
                    cat(cytoband_AH_assembly$assembly, sep = ", "), 
                    ". Please provide custom data for input assembly.",
                    call. = FALSE)
            
            cytoData <- NULL
        } else {
            
            # Load data from AHCytoBands
            if (!requireNamespace("AnnotationHub", 
                    quietly = TRUE)){
                
                warning("Please install `AnnotationHub` ",
                        "to plot an ideogram for a default",
                        " genome assembly.", call. = FALSE)
                cytoData <- NULL
                
            } else {
                
                assemblyName <- cytoband_AH_assembly[
                    cytoband_AH_assembly$assembly%in% assemblyName,]$assembly
                
                # Get name of AH object for assembly
                AH_id <- cytoband_AH_assembly[cytoband_AH_assembly$assembly == 
                                                assemblyName,]$AH
                
                # Check for internet
                if (has_internet()){

                    # Load AHCytoBands data
                    cytobands <- suppressMessages(AnnotationHub::query(
                        AnnotationHub::AnnotationHub(), "AHCytoBands"))
                    
                    # Grab data for assembly
                    cytoData <- as.data.frame(
                        suppressMessages(cytobands[[AH_id]]))
                } else {
                    # Try with localHub=TRUE for cache
                    errorFunction <- function(c){
                        return(NULL)
                    }
                    
                    cytobands <- 
                        tryCatch(
                            suppressMessages(
        AnnotationHub::query(AnnotationHub::AnnotationHub(localHub = TRUE),
                        "AHCytoBands")),
                    error = errorFunction
                        )
                    
                    if (!is.null(cytobands)){
                        cytoData <- as.data.frame(
                            suppressMessages(cytobands[[AH_id]]))
                    } else {
                        cytoData <- NULL
                    }
                    
                }
                
            }
            
        }

        return(cytoData)
    }

    ## Define a function to check that a chromosome name is in a TxDb
    checkChroms <- function(chrom, txdb) {
        if (is(txdb, "TxDb")) {
            tx_db <- txdb
        } else {
            tx_db <- eval(parse(text = paste0(as.name(txdb), "::", 
                                            as.name(txdb))))
        }

        txdbChroms <- GenomeInfoDb::seqlevels(tx_db)
        if (chrom %in% txdbChroms) {
            return(TRUE)
        } else {
            warning("'", chrom, "'",
                "not found in", txdb$packageName, ".",
                call. = FALSE
            )
            return(FALSE)
        }
    }

    ## Define a function that draws bands that fall within left curved regions
    curvedBands_left <- function(df, xCurve, yCurve, ymax) {
        start <- utils::type.convert(df["start"], as.is = TRUE)
        end <- utils::type.convert(df["end"], as.is = TRUE)
        col <- df["color"]

        if (end > max(xCurve)) {
            xpoints <- c(xCurve[which(xCurve >= start)], end, end)
            ypoints <- c(yCurve[which(xCurve >= start)], 0, ymax)
        } else {
            xpoints <- xCurve[which(xCurve >= start & xCurve <= end)]
            ypoints <- yCurve[which(xCurve >= start & xCurve <= end)]
        }

        if (length(xpoints) > 0 & length(ypoints) > 0) {
            curvedGrob <- polygonGrob(
                x = xpoints, y = ypoints,
                default.units = "native",
                gp = gpar(fill = col, col = NA)
            )

            assign("ideogram_grobs",
                addGrob(get("ideogram_grobs", envir = pgEnv),
                    child = curvedGrob
                ),
                envir = pgEnv
            )
        }
    }

    ## Define a function that draws bands that fall within right curved regions
    curvedBands_right <- function(df, xCurve, yCurve, ymax) {
        start <- utils::type.convert(df["start"], as.is = TRUE)
        end <- utils::type.convert(df["end"], as.is = TRUE)
        col <- df["color"]
        if (start < min(xCurve)) {
            xpoints <- c(xCurve[which(xCurve <= end)], start, start)
            ypoints <- c(yCurve[which(xCurve <= end)], ymax, 0)
        } else {
            xpoints <- xCurve[which(xCurve >= start & xCurve <= end)]
            ypoints <- yCurve[which(xCurve >= start & xCurve <= end)]
        }

        if (length(xpoints) > 0 & length(ypoints) > 0) {
            curvedGrob <- polygonGrob(
                x = xpoints, y = ypoints,
                default.units = "native",
                gp = gpar(fill = col, col = NA)
            )

            assign("ideogram_grobs",
                addGrob(get("ideogram_grobs", envir = pgEnv),
                    child = curvedGrob
                ),
                envir = pgEnv
            )
        }
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    ideoInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "ideoInternal"
    )
    
    ## Justification
    ideoInternal$just <- justConversion(just = ideoInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    ideogramPlot <- structure(list(
        chrom = ideoInternal$chrom,
        chromstart = 1, chromend = NULL,
        assembly = ideoInternal$assembly,
        colors = NULL,
        x = ideoInternal$x, y = ideoInternal$y,
        width = ideoInternal$width,
        height = ideoInternal$height,
        just = ideoInternal$just, grobs = NULL
    ),
    class = "ideogram"
    )
    attr(x = ideogramPlot, which = "plotted") <- ideoInternal$draw

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(ideogramPlot$chrom)) {
        stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
    }

    check_placement(object = ideogramPlot)
    errorcheck_plotIdeogram(orientation = ideoInternal$orientation)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    ideogramPlot$assembly <-
        parseAssembly(assembly = ideogramPlot$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    ideogramPlot <- defaultUnits(
        object = ideogramPlot,
        default.units = ideoInternal$default.units
    )

    # =========================================================================
    # GET APPROPRIATE BUILD DATA
    # =========================================================================

    if (is.null(ideoInternal$data)){
        data <- cytoAssembly(assembly = ideogramPlot$assembly)
    } else {
        data <- as.data.frame(ideoInternal$data)
        colnames(data) <- c("seqnames", "start", "end", "width", "strand",
                            "name", "gieStain")
    }
    
    ## TxDb data
    if (is(ideogramPlot$assembly$TxDb, "TxDb")){
        genome <- ideogramPlot$assembly$TxDb
    } else {
        if (!requireNamespace(ideogramPlot$assembly$TxDb, quietly = TRUE)){
            warning("`", ideogramPlot$assembly$TxDb, "` not available. ",
            "Please install to plot ideogram.", call. = FALSE)
            genome <- NULL
        } else {
            genome <- eval(parse(text = 
                        paste0(as.name(ideogramPlot$assembly$TxDb), "::",
                        as.name(ideogramPlot$assembly$TxDb))))
        }
    }
    
    chromLength <- 1
    if (!is.null(data) & !is.null(genome)) {
        
        chromCheck <- checkChroms(chrom = ideoInternal$chrom, 
                                txdb = ideogramPlot$assembly$TxDb)
        
        if (chromCheck == TRUE) {
            
            # =================================================================
            # ADD COLORS
            # =================================================================
            
            if (!is.null(ideoInternal$fill)){
                
                ## Check that length of fill vector is the same length as the
                ## of gieStain values in the data
                if (length(unique(data$gieStain)) != length(ideoInternal$fill)){
                    warning("`fill` color vector is not the same length as",
                            " number of gieStain values. Using default colors.")
                    colors <- cytobandColors
                } else {
                    
                    if (is.null(names(ideoInternal$fill))){
                        
                        ## No names, assign in order of levels  
                        data$gieStain <- as.factor(data$gieStain)
                        
                        colors <- data.frame("gieStain" = levels(data$gieStain),
                                             "color" = ideoInternal$fill)
                        
                    } else {
                        
                        ## Assign colors based on names
                        colors$gieStain <- 
                            data.frame("gieStain" = names(ideoInternal$fill),
                                       "color" = ideoInternal$fill)
                    }
                    
                } 
                
            } else {
                colors <- cytobandColors
            }
            
            data <- dplyr::left_join(data, colors, by = "gieStain")
            
            ## Add stains/colors to object
            objectCols <- unique(data[,c("gieStain", "color")])
            row.names(objectCols) <- NULL
            ideogramPlot$colors <- objectCols
            
            # =================================================================
            # SUBSET FOR CHROMOSOME
            # =================================================================

            data <- data[which(data[, 1] == ideogramPlot$chrom), ]
            data$seqnames <- as.character(data$seqnames)
            data$strand <- as.character(data$strand)
            data$name <- as.character(data$name)
            data$gieStain <- as.character(data$gieStain)
            chromLength <- GenomeInfoDb::seqlengths(genome)[[
            ideogramPlot$chrom]]
            ideogramPlot$chromend <- chromLength
        }
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(ideogramPlot$x) | is.null(ideogramPlot$y)) {
        height <- 0.075
        width <- 1

        scaleRatio <- width / height
        yscale <- chromLength / scaleRatio

        if (ideoInternal$orientation == "h") {
            vp <- viewport(
                height = unit(height, "snpc"),
                width = unit(width, "snpc"),
                x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                xscale = c(0, chromLength),
                yscale = c(0, yscale),
                just = "center",
                name = "ideogram1"
            )
        } else {
            height <- 1
            width <- 0.075
            vp <- viewport(
                height = unit(width, "snpc"),
                width = unit(height, "snpc"),
                x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                xscale = c(0, chromLength),
                yscale = c(0, yscale),
                angle = -90,
                just = "center",
                name = "ideogram1"
            )
        }

        if (ideoInternal$draw == TRUE) {
            grid.newpage()
        }
    } else {
        
        ## Get viewport name
        currentViewports <- current_viewports()
        vp_name <- paste0("ideogram", length(grep(
            pattern = "ideogram",
            x = currentViewports
        )) + 1)
        
        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = ideogramPlot)
        addViewport(vp_name)

        height <- convertHeight(page_coords$height,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
        width <- convertWidth(page_coords$width,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
        scaleRatio <- max(c(width, height)) / min(c(width, height))
        yscale <- chromLength / scaleRatio

        if (ideoInternal$orientation == "h") {
            vp <- viewport(
                height = page_coords$height,
                width = page_coords$width,
                x = page_coords$x, y = page_coords$y,
                xscale = c(0, chromLength),
                yscale = c(0, yscale),
                just = ideoInternal$just,
                name = vp_name
            )
        } else {
            ## Make viewport based on user inputs
            vpOG <- viewport(
                height = page_coords$height,
                width = page_coords$width,
                x = page_coords$x, y = page_coords$y,
                just = ideoInternal$just
            )

            ## Convert viewport to top right (top right of horizontal)
            vp_tr <- vp_topRight(viewport = vpOG)
            vp <- viewport(
                height = page_coords$width,
                width = page_coords$height,
                x = vp_tr[[1]], y = vp_tr[[2]],
                xscale = c(0, chromLength),
                yscale = c(0, yscale),
                angle = -90,
                just = c("left", "top"),
                name = vp_name
            )
        }
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS
    # =========================================================================

    assign("ideogram_grobs", gTree(vp = vp), envir = pgEnv)

    if (!is.null(data) & !is.null(genome)) {
        if (chromCheck == TRUE) {

            # =================================================================
            # CHROMOSOME GROBS
            # =================================================================

            ## Generate points along curves for the ends
            r <- vp$yscale[2] * 0.5
            leftAngles <- seq(pi / 2, 3 * pi / 2, pi / 500)
            rightAngles <- seq(3 * pi / 2, 5 * pi / 2, pi / 500)

            leftXpoints <- r + r * cos(leftAngles)
            leftYpoints <- r + r * sin(leftAngles)

            rightXpoints <- (chromLength - r) + r * cos(rightAngles)
            rightYpoints <- r + r * sin(rightAngles)


            if (nrow(data) > 1) {

                ## FIRST BAND ##
                firstBand <- data[which(data$start == 1), ]
                data <- subset(data, data$start != 1)

                if (firstBand$end > max(leftXpoints)) {
                    firstBand_Xpoints <- c(
                        leftXpoints, firstBand$end,
                        firstBand$end
                    )
                    firstBand_Ypoints <- c(leftYpoints, 0, vp$yscale[2])
                } else {
                    firstBand_Xpoints <- leftXpoints[which(leftXpoints <=
                        firstBand$end)]
                    firstBand_Ypoints <- leftYpoints[which(leftXpoints <=
                        firstBand$end)]
                }

                if (ideoInternal$showBands == TRUE) {
                    firstBand_grob <- polygonGrob(
                        x = firstBand_Xpoints,
                        y = firstBand_Ypoints,
                        default.units = "native",
                        gp = gpar(
                            fill = firstBand$color,
                            col = NA
                        )
                    )

                    assign("ideogram_grobs",
                        addGrob(get("ideogram_grobs", envir = pgEnv),
                            child = firstBand_grob
                        ),
                        envir = pgEnv
                    )
                }


                ## LAST BAND ##
                lastBand <- data[which(data$end == chromLength), ]
                data <- subset(data, data$end != chromLength)

                if (lastBand$start < min(rightXpoints)) {
                    lastBand_Xpoints <- c(
                        rightXpoints, lastBand$start,
                        lastBand$start
                    )
                    lastBand_Ypoints <- c(rightYpoints, vp$yscale[2], 0)
                } else {
                    lastBand_Xpoints <- rightXpoints[which(rightXpoints >=
                        lastBand$start)]
                    lastBand_Ypoints <- rightYpoints[which(rightXpoints >=
                        lastBand$start)]
                }
                if (ideoInternal$showBands == TRUE) {
                    lastBand_grob <- polygonGrob(
                        x = lastBand_Xpoints,
                        y = lastBand_Ypoints,
                        default.units = "native",
                        gp = gpar(
                            fill = lastBand$color,
                            col = NA
                        )
                    )

                    assign("ideogram_grobs",
                        addGrob(get("ideogram_grobs", envir = pgEnv),
                            child = lastBand_grob
                        ),
                        envir = pgEnv
                    )
                }



                if (ideoInternal$assembly %in% c("hg18", "hg19", "hg38")) {
                    ## CENTER BANDS ##
                    leftCent <- data[which(data$gieStain == "acen"), ][1, ]
                    leftCent_length <- leftCent$end - leftCent$start
                    rightCent <- data[which(data$gieStain == "acen"), ][2, ]
                    rightCent_length <- rightCent$end - rightCent$start
                    centerX <- leftCent$end
                    data <- subset(data, data$gieStain != "acen")

                    ## Generate points along curves for the centers
                    centerleftXpoints <- (centerX - r * 0.75) +
                        r * cos(rightAngles)
                    centerleftYpoints <- r + r * sin(rightAngles)
                    centerleftYpoints <- centerleftYpoints[which(
                        centerleftXpoints <= (rightCent$end - 0.5 *
                            rightCent_length)
                    )]
                    centerleftXpoints <- centerleftXpoints[which(
                        centerleftXpoints <= (rightCent$end - 0.5 *
                            rightCent_length)
                    )]

                    centerrightXpoints <- (centerX + r * 0.75) + r *
                        cos(leftAngles)
                    centerrightYpoints <- r + r * sin(leftAngles)
                    centerrightYpoints <- centerrightYpoints[which(
                        centerrightXpoints >= (leftCent$start + 0.5 *
                            leftCent_length)
                    )]
                    centerrightXpoints <- centerrightXpoints[which(
                        centerrightXpoints >= (leftCent$start + 0.5 *
                            leftCent_length)
                    )]

                    ## CENTER LEFT BAND ##
                    if (leftCent$start < min(centerleftXpoints)) {
                        leftCent_Xpoints <- c(
                            centerleftXpoints, leftCent$start,
                            leftCent$start
                        )
                        leftCent_Ypoints <- c(
                            centerleftYpoints,
                            vp$yscale[2], 0
                        )
                    } else {
                        leftCent_Xpoints <- centerleftXpoints[which(
                            centerleftXpoints >= leftCent$start
                        )]
                        leftCent_Ypoints <- centerleftYpoints[which(
                            centerleftXpoints >= leftCent$start
                        )]
                    }

                    if (ideoInternal$showBands == TRUE) {
                        leftCent_grob <- polygonGrob(
                            x = leftCent_Xpoints,
                            y = leftCent_Ypoints,
                            default.units = "native",
                            gp = gpar(fill = leftCent$color, col = NA)
                        )

                        assign("ideogram_grobs",
                            addGrob(get("ideogram_grobs", envir = pgEnv),
                                child = leftCent_grob
                            ),
                            envir = pgEnv
                        )
                    }


                    ## CENTER RIGHT BAND ##

                    if (rightCent$end > max(centerrightXpoints)) {
                        rightCent_Xpoints <- c(
                            centerrightXpoints, rightCent$end,
                            rightCent$end
                        )
                        rightCent_Ypoints <- c(
                            centerrightYpoints, 0,
                            vp$yscale[2]
                        )
                    } else {
                        rightCent_Xpoints <- centerrightXpoints[which(
                            centerrightXpoints <= rightCent$end
                        )]
                        rightCent_Ypoints <- centerrightYpoints[which(
                            centerrightXpoints <= rightCent$end
                        )]
                    }

                    if (ideoInternal$showBands == TRUE) {
                        rightCent_grob <- polygonGrob(
                            x = rightCent_Xpoints,
                            y = rightCent_Ypoints,
                            default.units = "native",
                            gp = gpar(
                                fill = rightCent$color,
                                col = NA
                            )
                        )

                        assign("ideogram_grobs",
                            addGrob(get("ideogram_grobs", envir = pgEnv),
                                child = rightCent_grob
                            ),
                            envir = pgEnv
                        )
                    }

                    ## GET ANY BANDS THAT FALL WITHIN CENTER CURVED REGIONS ##
                    if (ideoInternal$showBands == TRUE) {
                        inleftcurvedBands <- data[which(
                            data$end > min(centerleftXpoints) &
                                data$end <= centerX
                        ), ]
                        inrightcurvedBands <- data[which(
                            data$start < max(centerrightXpoints) &
                                data$start >= centerX
                        ), ]
                        if (nrow(inleftcurvedBands > 0)) {
                            invisible(apply(inleftcurvedBands,
                                1,
                                curvedBands_right,
                                xCurve = centerleftXpoints,
                                yCurve = centerleftYpoints,
                                ymax = vp$yscale[2]
                            ))
                        }
                        if (nrow(inrightcurvedBands > 0)) {
                            invisible(apply(inrightcurvedBands,
                                1,
                                curvedBands_left,
                                xCurve = centerrightXpoints,
                                yCurve = centerrightYpoints,
                                ymax = vp$yscale[2]
                            ))
                        }

                        ## REMAINING BANDS ##
                        data <- suppressMessages(dplyr::anti_join(
                            data, inleftcurvedBands
                        ))
                        data <- suppressMessages(dplyr::anti_join(
                            data, inrightcurvedBands
                        ))
                    }
                }


                if (ideoInternal$showBands == TRUE) {
                    ## GET ANY BANDS THAT FALL WITHIN OUTSIDE CURVED REGIONS ##
                    leftcurvedBands <- data[which(data$start <
                        max(leftXpoints)), ]
                    rightcurvedBands <- data[which(data$start >=
                        min(rightXpoints)), ]

                    if (nrow(leftcurvedBands > 0)) {
                        invisible(apply(leftcurvedBands,
                            1,
                            curvedBands_left,
                            xCurve = leftXpoints,
                            yCurve = leftYpoints,
                            ymax = vp$yscale[2]
                        ))
                    }
                    if (nrow(rightcurvedBands > 0)) {
                        invisible(apply(rightcurvedBands,
                            1,
                            curvedBands_right,
                            xCurve = rightXpoints,
                            yCurve = rightYpoints,
                            ymax = vp$yscale[2]
                        ))
                    }

                    ## REMAINING BANDS ##
                    data <- suppressMessages(dplyr::anti_join(
                        data,
                        leftcurvedBands
                    ))
                    data <- suppressMessages(dplyr::anti_join(
                        data,
                        rightcurvedBands
                    ))

                    rectBands <- rectGrob(
                        x = data$start, y = unit(0.5, "npc"),
                        width = data$width, height = unit(1, "npc"),
                        just = "left", default.units = "native",
                        gp = gpar(fill = data$color, col = NA)
                    )
                    assign("ideogram_grobs",
                        addGrob(get("ideogram_grobs", envir = pgEnv),
                            child = rectBands
                        ),
                        envir = pgEnv
                    )
                }
            }


            # =================================================================
            # OUTLINE GROBS
            # =================================================================

            if (ideoInternal$assembly %in% c("hg18", "hg19", "hg38")) {
                topIntersectY <- sqrt(r^2 - ((r * 0.75)^2)) + r
                bottomIntersectY <- -1 * sqrt(r^2 - ((r * 0.75)^2)) + r

                lbottomX <- centerleftXpoints[which(
                    centerleftXpoints <= centerX &
                        centerleftYpoints <= bottomIntersectY
                )]
                lbottomY <- centerleftYpoints[which(
                    centerleftXpoints <= centerX &
                        centerleftYpoints <= bottomIntersectY
                )]

                rbottomX <- centerrightXpoints[which(
                    centerrightXpoints >= centerX &
                        centerrightYpoints <= bottomIntersectY
                )]
                rbottomY <- centerrightYpoints[which(
                    centerrightXpoints >= centerX &
                        centerrightYpoints <= bottomIntersectY
                )]

                rtopX <- centerrightXpoints[which(
                    centerrightXpoints >= centerX &
                        centerrightYpoints >= topIntersectY
                )]
                rtopY <- centerrightYpoints[which(
                    centerrightXpoints >= centerX &
                        centerrightYpoints >= topIntersectY
                )]

                ltopX <- centerleftXpoints[which(
                    centerleftXpoints <= centerX &
                        centerleftYpoints >= topIntersectY
                )]
                ltopY <- centerleftYpoints[which(
                    centerleftXpoints <= centerX &
                        centerleftYpoints >= topIntersectY
                )]

                Xoutline <- c(
                    leftXpoints, lbottomX, centerX, rbottomX,
                    rightXpoints, rtopX, centerX, ltopX
                )
                Youtline <- c(
                    leftYpoints, lbottomY, bottomIntersectY, rbottomY,
                    rightYpoints, rtopY, topIntersectY, ltopY
                )
            } else {
                Xoutline <- c(leftXpoints, rightXpoints)
                Youtline <- c(leftYpoints, rightYpoints)
            }

            outlineGrob <- polygonGrob(
                x = Xoutline, y = Youtline,
                default.units = "native",
                gp = gpar(fill = NA, col = "#d0cfd4")
            )
            assign("ideogram_grobs",
                addGrob(get("ideogram_grobs", envir = pgEnv),
                    child = outlineGrob
                ),
                envir = pgEnv
            )
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (ideoInternal$draw == TRUE) {
        grid.draw(get("ideogram_grobs", envir = pgEnv))
    }
    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    ideogramPlot$grobs <- get("ideogram_grobs", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("ideogram[", vp$name, "]")
    invisible(ideogramPlot)
}
