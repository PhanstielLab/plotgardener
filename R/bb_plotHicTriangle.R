#' Plot a Hi-C interaction matrix in a triangular format
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
#' Default value is \code{assembly = "hg19"}.
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
#' data("bb_imrHicData")
#'
#' ## Create a page
#' bb_pageCreate(width = 4, height = 2.5, default.units = "inches")
#'
#' ## Plot and place triangle Hi-C plot
#' hicPlot <- bb_plotHicTriangle(data = bb_imrHicData, resolution = 10000,
#'                               zrange = c(0, 70),
#'                               chrom = "chr21",
#'                               chromstart = 28000000, chromend = 30300000,
#'                               x = 2, y = 0.5, width = 3, height = 1.5,
#'                               just = "top", default.units = "inches")
#'
#' ## Annotate x-axis genome label
#' bb_annoGenomeLabel(plot = hicPlot, scale = "Mb", x = 0.5, y = 2.03,
#'                    just = c("left", "top"))
#'
#' ## Annotate heatmap legend
#' bb_annoHeatmapLegend(plot = hicPlot, x = 3.5, y = 0.5,
#'                      width = 0.13, height = 1.2,
#'                      just = c("right", "top"))
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @details
#' A triangle Hi-C plot can be placed on a BentoBox coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' bb_plotHicTriangle(data, chrom,
#'                    chromstart = NULL, chromend = NULL,
#'                    x, y, width, height, just = c("left", "top"),
#'                    default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated triangle
#' Hi-C plot by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotHicTriangle(data, chrom,
#'                    chromstart = NULL, chromend = NULL)
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
                               assembly = "hg19",
                               palette = colorRampPalette(brewer.pal(n = 9, "YlGnBu")),
                               colorTrans = "linear", x = NULL, y = NULL,
                               width = NULL, height = NULL,
                               just = c("left", "top"),
                               default.units = "inches", draw = TRUE,
                               params = NULL, quiet = FALSE){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## For more accurate calculation of sqrt(2)
  two <- mpfr(2, 120)

  ## Define a function that catches errors for bb_plotTriangleHic
  errorcheck_bb_plotTriangleHic <- function(hic, hic_plot, norm, assembly){

    ###### hic/norm #####

    ## if it's a dataframe or datatable, it needs to be properly formatted
    if ("data.frame" %in% class(hic) && ncol(hic) != 3){

      stop("Invalid dataframe format.  Input a dataframe with 3 columns: chrA, chrB, counts.", call. = FALSE)

    }

    if (!"data.frame" %in% class(hic)){

      ## if it's a file path, it needs to be a .hic file
      if (file_ext(hic) != "hic"){

        stop("Invalid input. File must have a \".hic\" extension",
             call. = FALSE)

      }

      ## if it's a file path, it needs to exist
      if (!file.exists(hic)){

        stop(paste("File", hic, "does not exist."),
             call. = FALSE)

      }

      ## if it's a valid .hic file, it needs to have a valid norm parameter
      if (is.null(norm)){

        stop("If providing .hic file, please specify \'norm\'.",
             call. = FALSE)

      }

    }

    ##### chrom/chromstart/chromend #####


    ## Can't have only one NULL chromstart or chromend
    if ((is.null(hic_plot$chromstart)
         & !is.null(hic_plot$chromend))
        | (is.null(hic_plot$chromend)
           & !is.null(hic_plot$chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.",
           call. = FALSE)

    }

    ## Even though straw technically works without "chr" for hg19, will not accept for consistency purposes
    if (assembly == "hg19"){

      if (grepl("chr", hic_plot$chrom) == FALSE){

        stop(paste(paste0("'",hic_plot$chrom, "'"),
                   "is an invalid input for an hg19 chromsome. Please specify chromosome as",
                   paste0("'chr", hic_plot$chrom, "'.")),
             call. = FALSE)
      }

    }


    if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)){

      if (hic_plot$chromstart == hic_plot$chromend){
        stop("Genomic region is 0 bp long.", call. = FALSE)
      }

      ## Chromstart should be smaller than chromend
      if (hic_plot$chromstart > hic_plot$chromend){

        stop("\'chromstart\' should not be larger than \'chromend\'.",
             call. = FALSE)

      }

    }

    ##### zrange #####

    ## Ensure properly formatted zrange
    if (!is.null(hic_plot$zrange)){

      ## zrange needs to be a vector
      if (!is.vector(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be length 2
      if (length(hic_plot$zrange) != 2){

        stop("\'zrange\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of two numbers.",
             call. = FALSE)

      }

      ## second value should be larger than the first value
      if (hic_plot$zrange[1] >= hic_plot$zrange[2]){

        stop("\'zrange\' must be a vector of two numbers in which the 2nd value is larger than the 1st.", call. = FALSE)

      }

    }


    ##### height #####
    if (!is.null(hic_plot$height)){

      ## convert height to inches
      height <- convertHeight(hic_plot$height,
                              unitTo = "inches", valueOnly = TRUE)
      if (height < 0.05){
        stop("Height is too small for a valid triangle Hi-C plot.",
             call. = FALSE)

      }

    }

  }

  ## Define a function that subsets data
  subset_data <- function(hic, hic_plot){

    if (nrow(hic) > 0){
      hic <- hic[which(hic[,1] >= floor(hic_plot$chromstart/hic_plot$resolution)*hic_plot$resolution &
                         hic[,1] < hic_plot$chromend &
                         hic[,2] >= floor(hic_plot$chromstart/hic_plot$resolution)*hic_plot$resolution &
                         hic[,2] < hic_plot$chromend),]
    }


    return(hic)
  }

  ## Define a function that manually "clips" squares/triangles along edges
  manual_clip <- function(hic, hic_plot){

    clipLeft <- hic[which(hic[,1] < hic_plot$chromstart),]
    clipTop <- hic[which((hic[,2] + hic_plot$resolution) > hic_plot$chromend),]

    topLeft <- suppressMessages(dplyr::inner_join(clipLeft, clipTop))

    clipLeft <- suppressMessages(dplyr::anti_join(clipLeft, topLeft))
    clipTop <- suppressMessages(dplyr::anti_join(clipTop, topLeft))

    ############# Squares
    squares <- hic[which(hic[,2] > hic[,1]),]
    clipLeftsquares <- suppressMessages(dplyr::inner_join(squares, clipLeft))
    clipTopsquares <- suppressMessages(dplyr::inner_join(squares, clipTop))
    clippedSquares <- rbind(clipLeftsquares, clipTopsquares, topLeft)

    squares <- suppressMessages(dplyr::anti_join(squares, clippedSquares))


    clipLeftsquares$width <- hic_plot$resolution - (hic_plot$chromstart - clipLeftsquares$x)
    clipLeftsquares$x <- rep(hic_plot$chromstart, nrow(clipLeftsquares))

    #clipTopsquares$height <- hic_plot$resolution - (hic_plot$chromend - clipTopsquares$y)
    clipTopsquares$height <- hic_plot$chromend - clipTopsquares$y

    topLeft$width <-  hic_plot$resolution - (hic_plot$chromstart - topLeft$x)
    topLeft$x <- rep(hic_plot$chromstart, nrow(topLeft))
    #topLeft$height <- hic_plot$resolution - (hic_plot$chromend - topLeft$y)
    topLeft$height <- hic_plot$chromend - topLeft$y


    ############# Triangles
    triangles <- hic[which(hic[,2] == hic[,1]),]
    topRight <- suppressMessages(dplyr::inner_join(triangles, clipTop))
    bottomLeft <- suppressMessages(dplyr::inner_join(triangles, clipLeft))
    clippedTriangles <- rbind(topRight, bottomLeft)

    triangles <- suppressMessages(dplyr::anti_join(triangles, clippedTriangles))

    #topRight$height <- hic_plot$resolution - (hic_plot$chromend - topRight$y)
    topRight$height <- hic_plot$chromend - topRight$y
    topRight$width <- topRight$height

    bottomLeft$width <- hic_plot$resolution - (hic_plot$chromstart - bottomLeft$x)
    bottomLeft$height <- bottomLeft$width
    bottomLeft$x <- rep(hic_plot$chromstart, nrow(bottomLeft))
    bottomLeft$y <- rep(hic_plot$chromstart, nrow(bottomLeft))


    ## Recombine
    clippedHic <- rbind(squares, triangles, clipLeftsquares,
                        clipTopsquares, topLeft, topRight, bottomLeft)

    return(clippedHic)
  }

  ## Define a function that makes grobs for the triangle hic diagonal
  hic_diagonal <- function(hic){

    col <- hic[4]
    x <- as.numeric(hic[1])
    y <- as.numeric(hic[2])
    width <- as.numeric(hic[5])
    height <- as.numeric(hic[6])

    xleft = x
    xright = x + width
    ybottom = y
    ytop = y + height

    hic_triangle <- polygonGrob(x = c(xleft, xleft, xright),
                                y = c(ybottom, ytop, ytop),
                                gp = gpar(col = NA, fill = col),
                                default.units = "native")

    assign("hic_grobs2",
           addGrob(gTree = get("hic_grobs2", envir = bbEnv),
                   child = hic_triangle), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(resolution)) resolution <- NULL
  if(missing(palette)) palette <- NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(just)) just <- NULL
  if(missing(norm)) norm <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL
  if(missing(matrix)) matrix <- NULL
  if(missing(colorTrans)) colorTrans <- NULL
  if(missing(quiet)) quiet <- NULL

  ## Check if hic/chrom arguments are missing (could be in object)
  if(!hasArg(data)) data <- NULL
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_thicInternal <- structure(list(data = data, chrom = chrom,
                                    chromstart = chromstart,
                                    chromend = chromend,
                                    resolution = resolution,
                                    zrange = zrange, palette = palette,
                                    assembly = assembly, width = width,
                                    height = height, x = x,
                                    colorTrans = colorTrans,
                                    y = y, just = just, norm = norm,
                                    default.units = default.units,
                                    draw = draw, matrix = matrix,
                                    quiet = quiet),
                               class = "bb_thicInternal")

  bb_thicInternal <- parseParams(bb_params = params,
                                 object_params = bb_thicInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_thicInternal$resolution)) bb_thicInternal$resolution <- "auto"
  if(is.null(bb_thicInternal$palette)) bb_thicInternal$palette <- colorRampPalette(brewer.pal(n=9,"YlGnBu"))
  if(is.null(bb_thicInternal$assembly)) bb_thicInternal$assembly <- "hg19"
  if(is.null(bb_thicInternal$just)) bb_thicInternal$just <- c("left", "top")
  if(is.null(bb_thicInternal$norm)) bb_thicInternal$norm <- "KR"
  if(is.null(bb_thicInternal$default.units)) bb_thicInternal$default.units <- "inches"
  if(is.null(bb_thicInternal$draw)) bb_thicInternal$draw <- TRUE
  if(is.null(bb_thicInternal$matrix)) bb_thicInternal$matrix <- "observed"
  if(is.null(bb_thicInternal$colorTrans)) bb_thicInternal$colorTrans <- "linear"
  if(is.null(bb_thicInternal$quiet)) bb_thicInternal$quiet <- FALSE
  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  hic_plot <- structure(list(chrom = bb_thicInternal$chrom,
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
                             outsideVP = NULL, grobs = NULL),
                        class = "bb_hicTriangle")
  attr(x = hic_plot, which = "plotted") <- bb_thicInternal$draw

  # ======================================================================================================================================================================================
  # CHECK PLACEMENT/ARGUMENT ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_thicInternal$data)) stop("argument \"data\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_thicInternal$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
  check_placement(object = hic_plot)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  hic_plot$assembly <- parse_bbAssembly(assembly = hic_plot$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  hic_plot <- defaultUnits(object = hic_plot,
                           default.units = bb_thicInternal$default.units)

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotTriangleHic(hic = bb_thicInternal$data,
                                hic_plot = hic_plot,
                                norm = bb_thicInternal$norm,
                                assembly = hic_plot$assembly$Genome)

  # ======================================================================================================================================================================================
  # WHOLE CHROM INFORMATION
  # ======================================================================================================================================================================================

  if (is.null(hic_plot$chromstart) & is.null(hic_plot$chromend)){

    if (class(hic_plot$assembly$TxDb) == "TxDb"){
      txdbChecks <- TRUE
    } else {
      txdbChecks <- check_loadedPackage(package = hic_plot$assembly$TxDb,
                                        message = paste(paste0("`", hic_plot$assembly$TxDb,"`"),
                                                        "not loaded. Please install and load to plot full chromosome Hi-C map."))
    }

    scale <- c(0, 1)
    if (txdbChecks == TRUE){

      if (class(hic_plot$assembly$TxDb) == "TxDb"){
        tx_db <- hic_plot$assembly$TxDb
      } else {
        tx_db <- eval(parse(text = hic_plot$assembly$TxDb))
      }

      assembly_data <- GenomeInfoDb::seqlengths(tx_db)

      if (!hic_plot$chrom %in% names(assembly_data)){
        warning(paste("Chromosome",
                      paste0("'", hic_plot$chrom, "'"),
                      "not found in",
                      paste0("`", hic_plot$assembly$TxDb$packageName, "`"),
                      "and data for entire chromosome cannot be plotted."),
                call. = FALSE)
      } else {
        hic_plot$chromstart <- 1
        hic_plot$chromend <- assembly_data[[hic_plot$chrom]]
        hic_plot$altchromstart <- 1
        hic_plot$altchromend <- assembly_data[[hic_plot$chrom]]
        scale <- c(hic_plot$chromstart, hic_plot$chromend)
      }

    }

  } else {
    txdbChecks <- TRUE
    scale <- c(hic_plot$chromstart, hic_plot$chromend)
  }

  # ======================================================================================================================================================================================
  # ADJUST RESOLUTION
  # ======================================================================================================================================================================================

  if (bb_thicInternal$resolution == "auto"){
    hic_plot <- adjust_resolution(hic = bb_thicInternal$data,
                                  hic_plot = hic_plot)
  }

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  hic <- read_data(hic = bb_thicInternal$data, hic_plot = hic_plot,
                   norm = bb_thicInternal$norm, assembly = hic_plot$assembly,
                   type = bb_thicInternal$matrix,
                   quiet = bb_thicInternal$quiet)

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  hic <- subset_data(hic = hic, hic_plot = hic_plot)

  # ======================================================================================================================================================================================
  # SET ZRANGE AND SCALE DATA
  # ======================================================================================================================================================================================

  hic_plot <- set_zrange(hic = hic, hic_plot = hic_plot)
  hic$counts[hic$counts <= hic_plot$zrange[1]] <- hic_plot$zrange[1]
  hic$counts[hic$counts >= hic_plot$zrange[2]] <- hic_plot$zrange[2]

  # ======================================================================================================================================================================================
  # CONVERT NUMBERS TO COLORS
  # ======================================================================================================================================================================================

  ## if we don't have an appropriate zrange (even after setting it based on a null zrange), can't scale to colors
  if (!is.null(hic_plot$zrange) & length(unique(hic_plot$zrange)) == 2){

    if (grepl("log", bb_thicInternal$colorTrans) == TRUE){
      logBase <- as.numeric(gsub("log", "", bb_thicInternal$colorTrans))
      if (is.na(logBase)){
        logBase <- exp(1)
      }

      ## Won't scale to log if negative values
      if (any(hic$counts < 0)){
        stop("Negative values in Hi-C data. Cannot scale colors on a log scale. Please
             set `colorTrans = 'linear'`.", call. = FALSE)
      }


      hic$counts <- log(hic$counts, base = logBase)
      hic$color <- bb_maptocolors(hic$counts,
                                  col = bb_thicInternal$palette,
                                  num = 100,
                                  range = log(hic_plot$zrange, logBase))
    } else {
      hic$color <- bb_maptocolors(hic$counts,
                                  col = bb_thicInternal$palette,
                                  num = 100,
                                  range = hic_plot$zrange)
    }

    hic_plot$color_palette <- bb_thicInternal$palette

    }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_hicTriangle",
                    length(grep(pattern = "bb_hicTriangle",
                                x = currentViewports)) + 1)

  if (is.null(hic_plot$x) & is.null(hic_plot$y)){

    inside_vp <- viewport(height = unit(1, "npc"), width = unit(0.5, "npc"),
                          x = unit(0, "npc"), y = unit(0, "npc"),
                          xscale = scale,
                          yscale = scale,
                          just = c("left", "bottom"),
                          name = paste0(vp_name, "_inside"),
                          angle = -45)

    outside_vp <- viewport(height = unit(0.75, "snpc"),
                           width = unit(1.5, "snpc"),
                           x = unit(0.125, "npc"),
                           y = unit(0.25, "npc"),
                           xscale = scale,
                           clip = "on",
                           just = c("left", "bottom"),
                           name = paste0(vp_name, "_outside"))


    if (bb_thicInternal$draw == TRUE){

      vp_name <- "bb_hicTriangle1"
      inside_vp$name <- "bb_hicTriangle1_inside"
      outside_vp$name <- "bb_hicTriangle1_outside"
      grid.newpage()

    }

  } else {

    ## Get sides of viewport based on input width
    vp_side <- (convertWidth(hic_plot$width,
                             unitTo = get("page_units", envir = bbEnv),
                             valueOnly = TRUE))/sqrt(two)

    ## Convert coordinates into same units as page for outside vp
    page_coords <- convert_page(object = hic_plot)

    ## Get bottom left point of triangle (hence bottom left of actual viewport) based on just
    bottom_coords <- vp_bottomLeft(viewport(x = page_coords$x,
                                            y = page_coords$y,
                                            width = page_coords$width,
                                            height = page_coords$height,
                                            just = hic_plot$just))

    inside_vp <- viewport(height = unit(vp_side,
                                        get("page_units", envir = bbEnv)),
                          width = unit(vp_side,
                                       get("page_units", envir = bbEnv)),
                          x = unit(0, "npc"),
                          y = unit(0, "npc"),
                          xscale = scale,
                          yscale = scale,
                          just = c("left", "bottom"),
                          name = paste0(vp_name, "_inside"),
                          angle = -45)

    outside_vp <- viewport(height = page_coords$height,
                           width = page_coords$width,
                           x = unit(bottom_coords[[1]],
                                    get("page_units", envir = bbEnv)),
                           y = unit(bottom_coords[[2]],
                                    get("page_units", envir = bbEnv)),
                           xscale = scale,
                           clip = "on",
                           just = c("left", "bottom"),
                           name = paste0(vp_name, "_outside"))

    add_bbViewport(paste0(vp_name, "_outside"))
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  hic_plot$outsideVP <- outside_vp
  assign("hic_grobs2", gTree(vp = inside_vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)){

    if (nrow(hic) > 0){

      hic$width <- hic_plot$resolution
      hic$height <- hic_plot$resolution

      ## Manually "clip" the grobs that fall out of the desired chromstart to chromend region
      hic <- manual_clip(hic = hic, hic_plot = hic_plot)
      hic <- hic[order(as.numeric(rownames(hic))),]

      ## Separate into squares for upper region and triangle shapes for the diagonal
      squares <- hic[which(hic[,2] > hic[,1]),]
      triangles <- hic[which(hic[,2] == hic[,1]),]

      if (nrow(squares) > 0){

        ## Make square grobs and add to grob gTree
        hic_squares <- rectGrob(x = squares$x,
                                y = squares$y,
                                just = c("left", "bottom"),
                                width = squares$width,
                                height = squares$height,
                                gp = gpar(col = NA, fill = squares$color),
                                default.units = "native")
        assign("hic_grobs2",
               addGrob(gTree = get("hic_grobs2", envir = bbEnv),
                       child = hic_squares), envir = bbEnv)
      }

      if (nrow(triangles) > 0){
        ## Make triangle grobs and add to grob gTree
        invisible(apply(triangles, 1, hic_diagonal))

      }

      if (nrow(squares) == 0 & nrow(triangles) == 0){

        if (txdbChecks == TRUE){
          warning("Warning: no data found in region.  Suggestions: check chromosome, check region.", call. = FALSE)
        }

    }


    } else {

      warning("Warning: no data found in region.  Suggestions: check chromosome, check region.", call. = FALSE)

    }


  }


  # ======================================================================================================================================================================================
  # IF DRAW == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_thicInternal$draw == TRUE){

    pushViewport(outside_vp)
    grid.draw(get("hic_grobs2", envir = bbEnv))
    upViewport()
  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  hic_plot$grobs <- get("hic_grobs2", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_hicTriangle[", vp_name, "]"))
  invisible(hic_plot)

}
