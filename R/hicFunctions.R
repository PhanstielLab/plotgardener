## Define a function to adjust/detect resolution based on .hic file/dataframe
# @param hic hic data argument
# @param hic_plot hic plot object
adjust_resolution <- function(hic, hic_plot) {
    if (!("data.frame" %in% class(hic))) {
        if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)) {
            fileResolutions <- strawr::readHicBpResolutions(hic)

            ## Get range of data and try to pick a resolution
            dataRange <- hic_plot$chromend - hic_plot$chromstart
            if (dataRange >= 150000000) {
                bestRes <- max(fileResolutions)
            } else if (dataRange >= 75000000 & dataRange < 150000000) {
                bestRes <- 250000
                bestRes <- fileResolutions[which(
                    abs(fileResolutions - bestRes) == min(
                        abs(fileResolutions - bestRes)
                    )
                )]
            } else if (dataRange >= 35000000 & dataRange < 75000000) {
                bestRes <- 100000
                bestRes <- fileResolutions[which(
                    abs(fileResolutions - bestRes) == min(
                        abs(fileResolutions - bestRes)
                    )
                )]
            } else if (dataRange >= 20000000 & dataRange < 35000000) {
                bestRes <- 50000
                bestRes <- fileResolutions[which(
                    abs(fileResolutions - bestRes) == min(
                        abs(fileResolutions - bestRes)
                    )
                )]
            } else if (dataRange >= 5000000 & dataRange < 20000000) {
                bestRes <- 25000
                bestRes <- fileResolutions[which(
                    abs(fileResolutions - bestRes) == min(
                        abs(fileResolutions - bestRes)
                    )
                )]
            } else if (dataRange >= 3000000 & dataRange < 5000000) {
                bestRes <- 10000
                bestRes <- fileResolutions[which(
                    abs(fileResolutions - bestRes) == min(
                        abs(fileResolutions - bestRes)
                    )
                )]
            } else {
                bestRes <- 5000
                bestRes <- fileResolutions[which(
                    abs(fileResolutions - bestRes) == min(
                        abs(fileResolutions - bestRes)
                    )
                )]
            }

            hic_plot$resolution <- as.integer(bestRes)
        }
    } else {

        ## Try to detect resolution from data
        offDiag <- hic[which(hic[, 1] != hic[, 2]), ]
        bpDiffs <- abs(offDiag[, 2] - offDiag[, 1])
        predRes <- min(bpDiffs)

        hic_plot$resolution <- as.integer(predRes)
    }

    return(hic_plot)
}

## Define a function that detects bin limit for plotting Hi-C data
# @param hic hic data argument
# @param hic_plot hic plot object
hic_limit <- function(hic, hic_plot){
    ## Calculate bin number
    dataRange <- hic_plot$chromend - hic_plot$chromstart
    binNumber <- dataRange/hic_plot$resolution
    
    if (binNumber > 1000){
        
        if (!("data.frame" %in% class(hic))) {
            
            ## Overwrite manual resolution for Hi-C file
            hic_plotNew <- adjust_resolution(hic = hic, hic_plot = hic_plot)
            newRes <- hic_plotNew$resolution
            hic_plot$resolution <- newRes
            warning("Attempting to plot too many Hi-C pixels. Adjusting to ",
                    "a resolution of ", newRes, " BP.", call. = FALSE)
        } else {
            warning(hic_plot$resolution, " BP resolution detected in ",
            "dataframe. Attempting to plot too many Hi-C pixels. Please ",
            "read in data at a lower resolution before attempting to plot.", 
            call. = FALSE)
            hic_plot$resolution <- NA
        }
    
    }
    
    return(hic_plot)
}

## Define a function that reads in hic data for plotHic functions
# @param hic hic data argument
# @param hic_plot hic plot object
# @param norm normalization factor
# @param assembly genome assembly
# @param type matrix type
# @param quiet message quiet parameter
read_data <- function(hic, hic_plot, norm, assembly, type, quiet) {

    ## if .hic file, read in with bb_rhic
    if (!("data.frame" %in% class(hic))) {
        if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend) &
            !is.na(hic_plot$resolution)) {
            readchromstart <- hic_plot$chromstart - hic_plot$resolution
            readchromend <- hic_plot$chromend + hic_plot$resolution
            readaltchromstart <- hic_plot$altchromstart - hic_plot$resolution
            readaltchromend <- hic_plot$altchromend + hic_plot$resolution
            hic <- suppressWarnings(bb_readHic(
                file = hic, chrom = hic_plot$chrom,
                chromstart = readchromstart,
                chromend = readchromend,
                assembly = assembly,
                resolution = hic_plot$resolution,
                zrange = hic_plot$zrange,
                norm = norm,
                altchrom = hic_plot$altchrom,
                altchromstart = readaltchromstart,
                altchromend = readaltchromend,
                matrix = type
            ))
        } else {
            hic <- data.frame(matrix(nrow = 0, ncol = 3))
        }
    } else {
    
        if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend) &
            !is.na(hic_plot$resolution)) {
            if (!quiet) {
                message(
                    "Read in dataframe.  Assuming \'chrom\' in column1 ",
                    "and \'altchrom\' in column2. ",
                    hic_plot$resolution, " BP resolution detected."
                )
            }
            ## bb_plotHicRectangle specific warning for missing data
            if (hic_plot$chromstart < min(hic[, 1]) |
                hic_plot$chromend > max(hic[, 1])) {
                warning("`bb_plotHicRectangle` requires additional data to",
                " plot a rectangular plot. Data is missing from input ",
                "dataframe for region and plot will be a trapezoid. To avoid ",
                "this missing data, call `bb_plotHicRectangle` with full .hic",
                " file.", call. = FALSE
                )
            }
        } else {
            hic <- data.frame(matrix(nrow = 0, ncol = 3))
        }
    }

    ## Rename columns for later processing
    colnames(hic) <- c("x", "y", "counts")
    hic <- na.omit(hic)

    return(hic)
}

## Define a function that sets the Hi-C zrange
# @param hic hic data argument
# @param hic_plot hic plot object
set_zrange <- function(hic, hic_plot) {

    ## no zrange, only one value
    if (is.null(hic_plot$zrange) & length(unique(hic$counts)) == 1) {
        zrange <- c(unique(hic$counts), unique(hic$counts))
        hic_plot$zrange <- zrange
    }

    ## no zrange, multiple values
    if (is.null(hic_plot$zrange) & length(unique(hic$counts)) > 1) {
        if (grepl("log", hic_plot$colorTrans) == TRUE) {
            zrange <- c(0.0001, max(hic$counts))
        } else {
            zrange <- c(0, max(hic$counts))
        }

        hic_plot$zrange <- zrange
    }

    return(hic_plot)
}

## Define a function that parses an inherited half of a Hi-C plot
# @param hic hic plot object
inherit_half <- function(hic) {
    if (is(hic, "bb_hicSquare")) {
        half <- hic$half
    } else if (is(hic, "bb_hicTriangle") |
        is(hic, "bb_hicRectangle")) {
        half <- "top"
    }

    return(half)
}