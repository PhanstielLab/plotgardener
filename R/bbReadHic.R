#' Read a .hic file and return Hi-C data as a dataframe
#' 
#' @usage bbReadHic(
#'     file,
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     altchrom = NULL,
#'     altchromstart = NULL,
#'     altchromend = NULL,
#'     assembly = "hg38",
#'     resolution = "auto",
#'     res_scale = "BP",
#'     zrange = NULL,
#'     norm = "KR",
#'     matrix = "observed",
#'     params = NULL,
#'     quiet = FALSE
#' )
#'
#' @param file A character value specifying the path to the .hic file.
#' @param chrom Chromosome of data, as a string.
#' @param chromstart Integer start position on chromosome.
#' @param chromend Integer end position on chromosome.
#' @param altchrom Alternate chromosome for interchromosomal data,
#' as a string.
#' @param altchromstart Alternate chromosome integer start position
#' for interchromosomal data.
#' @param altchromend Alternate chromosome integer end position
#' for interchromosomal data.
#' @param assembly Default genome assembly as a string or a
#' \link[BentoBox]{bbAssembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param resolution A numeric specifying the width of each pixel.
#' "auto" will attempt to choose a resolution in basepairs based on
#' the size of the region.
#' @param res_scale A character value specifying the resolution scale.
#' Default value is \code{res_scale = "BP"}. Options are:
#' \itemize{
#' \item{\code{"BP"}: }{Base pairs.}
#' \item{\code{"FRAG"}: }{Fragments.}
#' }
#' @param zrange A numeric vector of length 2 specifying the range of
#' interaction scores, where extreme values will be set to the max or min.
#' @param norm Character value specifying hic data normalization method.
#' This value must be found in the .hic file.
#' Default value is \code{norm = "KR"}.
#' @param matrix Character value indicating the type of matrix to output.
#' Default value is \code{matrix = "observed"}. Options are:
#' \itemize{
#' \item{\code{"observed"}: }{Observed counts.}
#' \item{\code{"oe"}: }{Observed/expected counts.}
#' \item{\code{"log2oe"}: }{Log2 transformed observed/expected counts.}
#' }
#' @param params An optional \link[BentoBox]{bbParams} object
#' containing relevant function parameters.
#' @param quiet A logical indicating whether or not to print messages.
#'
#' @return Returns a 3-column dataframe in sparse upper triangular
#' format with the following columns: \code{chrom}, \code{altchrom},
#' \code{counts}.
#' 
#' @examples 
#' hicFile <- system.file("extdata/test_chr22.hic", package="BentoBoxData")
#' 
#' ## Read in data for all chr22 file at 2.5Mb bp resolution
#' hicData <- bbReadHic(file = hicFile, chrom = "chr22",
#'                     assembly = "hg19",
#'                     resolution = 2500000) 
#'                         
#' ## Read in region `chr22:20000000-47500000` at 100 Kb resolution
#' hicData10Kb <- bbReadHic(file = hicFile, chrom = "chr22",
#'                         chromstart = 20000000, chromend = 47500000,
#'                         assembly = "hg19",
#'                         resolution = 100000)                     
#'
#' @seealso \link[strawr]{straw}
#'
#' @export
bbReadHic <- function(file, chrom, chromstart = NULL, chromend = NULL,
                    altchrom = NULL, altchromstart = NULL,
                    altchromend = NULL, assembly = "hg38",
                    resolution = "auto", res_scale = "BP",
                    zrange = NULL, norm = "KR", matrix = "observed",
                    params = NULL, quiet = FALSE) {


    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for bb_rhic
    errorcheck_bb_rhic <- function(hic, chrom, chromstart, chromend,
                                zrange, altchrom, altchromstart,
                                altchromend, norm, res_scale, assembly) {

        ## hic input needs to be a path to a .hic file
        if (!is(hic, "character")) {
            stop("Invalid input. Input needs to be a path to a .hic file.",
                call. = FALSE
            )
        }

        if ((file_ext(hic) != "hic")) {
            stop("Invalid input. File must have a \".hic\" extension",
                call. = FALSE
            )
        }

        if (!file.exists(hic)) {
            stop("File", hic, "does not exist.", call. = FALSE)
        }

        ## Not supporting chrM
        if (chrom == "chrM") {
            stop("chrM not supported.", call. = FALSE)
        }

        ## Even though straw technically works without "chr" for hg19,
        ## will not accept for consistency purposes
        if (assembly == "hg19") {
            if (grepl("chr", chrom) == FALSE) {
                stop("'", chrom, "'",
                    "is an invalid input for an hg19 ",
                    "chromsome. Please specify chromosome as",
                    "'chr", chrom, "'.",
                    call. = FALSE
                )
            }
        }

        ## Genomic region
        bb_regionErrors(chromstart = chromstart,
                    chromend = chromend)
        
        if (!is.null(altchrom)) {
            if (altchrom == "chrM") {
                stop("chrM not supported.", call. = FALSE)
            }

            ## Even though straw technically works without "chr" for hg19,
            ## will not accept for consistency purposes
            if (assembly == "hg19") {
                if (grepl("chr", altchrom) == FALSE) {
                    stop("'", altchrom, "'",
                        "is an invalid input for an hg19 chromsome. ",
                        "Please specify chromosome as",
                        "'chr", altchrom, "'.",
                        call. = FALSE
                    )
                }
            }

            ## Can't specify altchrom without a chrom
            if (is.null(chrom)) {
                stop("Specified \'altchrom\', but did not give \'chrom\'.",
                    call. = FALSE
                )
            }
            
            bb_regionErrors(chromstart = altchromstart,
                        chromend = altchromend)
            
            ## If giving same chrom and altchrom, need to specify
            ## chromstart/chromend and altchromstart/altchromend
            
            if (chrom == altchrom) {
                if (is.null(chromstart) |
                    is.null(chromend) |
                    is.null(altchromstart) |
                    is.null(altchromend)) {
                    stop("If giving the same \'chrom\' and \'altchrom\', ",
                        "please specify \'chromstart\', \'chromend\', ",
                        "\'altchromstart\', and \'altchromend\'. ",
                        "If trying to get all interactions between one ",
                        "chromosome, just specify \'chrom\'.", call. = FALSE)
                }
            }
        }

        ## zrange errors
        bb_rangeErrors(range = zrange)

        ## Check for valid "res_scale" parameter
        if (!(res_scale %in% c("BP", "FRAG"))) {
            stop("Invalid \'res_scale\'.  Options are \'BP\' and \'FRAG\'.",
                call. = FALSE
            )
        }
    }

    ## Define a function that determines a best resolution for size of region
    auto_resolution <- function(file, chromstart, chromend) {
        fileResolutions <- strawr::readHicBpResolutions(file)
        if (is.null(chromstart) & is.null(chromend)) {
            autoRes <- max(fileResolutions)
        } else {
            dataRange <- chromend - chromstart
            if (dataRange >= 150000000) {
                autoRes <- max(fileResolutions)
            } else if (dataRange >= 75000000 & dataRange < 150000000) {
                autoRes <- 250000
                autoRes <- fileResolutions[which(
                    abs(fileResolutions - autoRes) == min(
                        abs(fileResolutions - autoRes)
                    )
                )]
            } else if (dataRange >= 35000000 & dataRange < 75000000) {
                autoRes <- 100000
                autoRes <- fileResolutions[which(
                    abs(fileResolutions - autoRes) == min(
                        abs(fileResolutions - autoRes)
                    )
                )]
            } else if (dataRange >= 20000000 & dataRange < 35000000) {
                autoRes <- 50000
                autoRes <- fileResolutions[which(
                    abs(fileResolutions - autoRes) == min(
                        abs(fileResolutions - autoRes)
                    )
                )]
            } else if (dataRange >= 5000000 & dataRange < 20000000) {
                autoRes <- 25000
                autoRes <- fileResolutions[which(
                    abs(fileResolutions - autoRes) == min(
                        abs(fileResolutions - autoRes)
                    )
                )]
            } else if (dataRange >= 3000000 & dataRange < 5000000) {
                autoRes <- 10000
                autoRes <- fileResolutions[which(
                    abs(fileResolutions - autoRes) == min(
                        abs(fileResolutions - autoRes)
                    )
                )]
            } else {
                autoRes <- 5000
                autoRes <- fileResolutions[which(
                    abs(fileResolutions - autoRes) == min(
                        abs(fileResolutions - autoRes)
                    )
                )]
            }
        }

        return(as.integer(autoRes))
    }

    ## Define a function to parse chromsome/region for Straw
    parse_region <- function(chrom, chromstart, chromend, assembly) {
        if (assembly == "hg19") {
            strawChrom <- gsub("chr", "", chrom)
        } else {
            strawChrom <- chrom
        }

        if (is.null(chromstart) & is.null(chromend)) {
            regionStraw <- strawChrom
        } else {

            ## Keep chromstart and chromend without scientific notation
            ## for processing with Straw

            chromstart <- format(chromstart, scientific = FALSE)
            chromend <- format(chromend, scientific = FALSE)
            regionStraw <- paste(strawChrom, chromstart, sep = ":")
            regionStraw <- paste(regionStraw, chromend, sep = "-")
        }

        return(regionStraw)
    }

    ## Define a function to reorder chromsomes to put "chrom" input in col1
    orderChroms <- function(hic, chrom, altchrom, assembly) {
        
        chrom <- gsub("chr", "", chrom)
        altchrom <- gsub("chr", "", altchrom)
        
        if (!"X" %in% chrom & !"Y" %in% chrom) {
            chrom <- utils::type.convert(chrom, as.is = TRUE)
        }

        if (!"X" %in% altchrom & !"Y" %in% altchrom) {
            altchrom <- utils::type.convert(altchrom, as.is = TRUE)
        }

        ## CASE 1: two numbers
        if (all(is(chrom, "numeric"), is(altchrom, "numeric"))) {
            if (chrom > altchrom) {
                hic <- hic[, c("y", "x", "counts")]
            }
        } else if (any(is(chrom, "numeric"), is(altchrom, "numeric"))) {
            ## CASE 2: number and X/Y
            if (is(altchrom, "numeric")) {
                hic <- hic[, c("y", "x", "counts")]
            }
        } else {
            ## CASE 3: X and Y
            if ("Y" %in% chrom) {
                hic <- hic[, c("y", "x", "counts")]
            }
        }

        return(hic)
    }

    ## Define a function to scale data with zrange
    scale_data <- function(upper, zrange) {
        if (!is.null(zrange)) {
            upper$counts[upper$counts <= zrange[1]] <- zrange[1]
            upper$counts[upper$counts >= zrange[2]] <- zrange[2]
        } else {

            # if null, zrange will be set to (0, max(data))
            upper$counts[upper$counts <= 0] <- 0
        }

        return(upper)
    }

    ## Define a function to rename columns
    rename_columns <- function(upper, chrom, altchrom) {
        if (is.null(altchrom)) {
            colnames(upper) <- c(paste0(chrom, "_A"), 
                            paste0(chrom, "_B"), "counts")
        } else {
            colnames(upper) <- c(chrom, altchrom, "counts")
        }

        return(upper)
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_rhic <- parseParams(params = params, 
                        defaultArgs = formals(eval(match.call()[[1]])),
                        declaredArgs = lapply(match.call()[-1], 
                                            eval.parent, n = 2),
                        class = "bb_rhic")
    
    if (is.null(bb_rhic$file)) stop("argument \"file\" is missing, ",
                                    "with no default.", call. = FALSE)
    if (is.null(bb_rhic$chrom)) stop("argument \"chrom\" is missing, ",
                                    "with no default.", call. = FALSE)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    bb_rhic$assembly <- parse_bbAssembly(assembly = bb_rhic$assembly)

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    errorcheck_bb_rhic(
        hic = bb_rhic$file, chrom = bb_rhic$chrom,
        chromstart = bb_rhic$chromstart,
        chromend = bb_rhic$chromend, zrange = bb_rhic$zrange,
        altchrom = bb_rhic$altchrom,
        altchromstart = bb_rhic$altchromstart,
        altchromend = bb_rhic$altchromend, norm = bb_rhic$norm,
        res_scale = bb_rhic$res_scale,
        assembly = bb_rhic$assembly$Genome
    )

    # =========================================================================
    # SET PARAMETERS
    # =========================================================================

    parse_chromstart <- bb_rhic$chromstart
    parse_chromend <- bb_rhic$chromend
    parse_altchromstart <- bb_rhic$altchromstart
    parse_altchromend <- bb_rhic$altchromend


    ## For off diagonal plotting, grabbing whole symmetric region
    if (!is.null(bb_rhic$altchrom)) {
        if (bb_rhic$chrom == bb_rhic$altchrom) {
            parse_chromstart <- min(bb_rhic$chromstart, bb_rhic$altchromstart)
            parse_chromend <- max(bb_rhic$chromend, bb_rhic$altchromend)
        }
    }

    # =========================================================================
    # PARSE REGIONS
    # =========================================================================

    chromRegion <- parse_region(
        chrom = bb_rhic$chrom,
        chromstart = parse_chromstart,
        chromend = parse_chromend,
        assembly = bb_rhic$assembly$Genome
    )

    if (is.null(bb_rhic$altchrom)) {
        altchromRegion <- chromRegion
    } else {
        if (bb_rhic$chrom == bb_rhic$altchrom) {
            altchromRegion <- chromRegion
        } else {
            altchromRegion <- parse_region(
                chrom = bb_rhic$altchrom,
                chromstart = parse_altchromstart,
                chromend = parse_altchromend,
                assembly = bb_rhic$assembly$Genome
            )
        }
    }

    # =========================================================================
    # ADJUST RESOLUTION
    # =========================================================================

    if (bb_rhic$resolution == "auto") {
        bb_rhic$resolution <- auto_resolution(
            file = bb_rhic$file,
            chromstart = bb_rhic$chromstart,
            chromend = bb_rhic$chromend
        )
        bb_rhic$res_scale <- "BP"
    }

    # =========================================================================
    # EXTRACT SPARSE UPPER TRIANGULAR USING STRAW
    # =========================================================================

    errorFunction <- function(c) {
        upper <- data.frame(matrix(nrow = 0, ncol = 3))
        colnames(upper) <- c("x", "y", "counts")
        return(upper)
    }

    log <- FALSE
    if (bb_rhic$matrix == "logoe") {
        bb_rhic$matrix <- "oe"
        log <- TRUE
    }
    
    print(bb_rhic$norm)
    print(bb_rhic$file)
    print(toString(chromRegion))
    print(toString(altchromRegion))
    print(bb_rhic$res_scale)
    print(bb_rhic$resolution)
    print(bb_rhic$matrix)
    
    upper <-
        tryCatch(strawr::straw(
            norm = bb_rhic$norm,
            fname = bb_rhic$file,
            chr1loc = toString(chromRegion),
            chr2loc = toString(altchromRegion),
            unit = bb_rhic$res_scale,
            binsize = bb_rhic$resolution,
            matrix = bb_rhic$matrix
        ),
        error = errorFunction
        )

    if (log == TRUE) {
        upper[, "counts"] <- log2(upper[, "counts"])
    }


    # =========================================================================
    # REORDER COLUMNS BASED ON CHROM/ALTCHROM INPUT
    # =========================================================================

    if (!is.null(bb_rhic$altchrom)) {
        upper <- orderChroms(
            hic = upper, chrom = bb_rhic$chrom,
            altchrom = bb_rhic$altchrom,
            assembly = bb_rhic$assembly$Genome
        )
        colnames(upper) <- c("x", "y", "counts")
    }

    # =========================================================================
    # SCALE DATA WITH ZRANGE
    # =========================================================================

    scaled_data <- scale_data(upper = upper, zrange = bb_rhic$zrange)

    # =========================================================================
    # FORMAT DATA IN PROPER ORDER AND WITH LABELS
    # =========================================================================

    renamed_data <- rename_columns(
        upper = scaled_data,
        chrom = bb_rhic$chrom,
        altchrom = bb_rhic$altchrom
    )

    # =========================================================================
    # REMOVE NAN VALUES
    # =========================================================================

    renamed_data <- na.omit(renamed_data)

    # =========================================================================
    # RETURN DATAFRAME
    # =========================================================================
    if (nrow(renamed_data) == 0) {
        warning("No data found in region. Suggestions: check that ",
                "chromosome names match genome assembly with ",
                "`strawr::readHicChroms()`; ",
                "check region; check available resolutions ",
                "with `strawr::readHicBpResolutions()`.", call. = FALSE)
    } else {
        if (!bb_rhic$quiet) {
            message(
                "Read in hic file with ",
                bb_rhic$norm, " normalization at ",
                bb_rhic$resolution, " ", bb_rhic$res_scale,
                " resolution."
            )
        }
    }
    return(renamed_data)
}
