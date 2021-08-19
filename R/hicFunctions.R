## Define a function to adjust/detect resolution based on .hic file/dataframe
# @param hic hic data argument
# @param hicPlot hic plot object
adjust_resolution <- function(hic, hicPlot) {
    if (!("data.frame" %in% class(hic))) {
        if (!is.null(hicPlot$chromstart) & !is.null(hicPlot$chromend)) {
            fileResolutions <- strawr::readHicBpResolutions(hic)

            ## Get range of data and try to pick a resolution
            dataRange <- hicPlot$chromend - hicPlot$chromstart
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

            hicPlot$resolution <- as.integer(bestRes)
        }
    } else {

        ## Try to detect resolution from data
        offDiag <- hic[which(hic[, 1] != hic[, 2]), ]
        bpDiffs <- abs(offDiag[, 2] - offDiag[, 1])
        predRes <- min(bpDiffs)

        hicPlot$resolution <- as.integer(predRes)
    }

    return(hicPlot)
}

## Define a function that detects bin limit for plotting Hi-C data
# @param hic hic data argument
# @param hicPlot hic plot object
hic_limit <- function(hic, hicPlot){
    ## Calculate bin number
    dataRange <- hicPlot$chromend - hicPlot$chromstart
    binNumber <- dataRange/hicPlot$resolution
    
    if (binNumber > 1000){
        
        if (!("data.frame" %in% class(hic))) {
            
            ## Overwrite manual resolution for Hi-C file
            hicPlotNew <- adjust_resolution(hic = hic, hicPlot = hicPlot)
            newRes <- hicPlotNew$resolution
            hicPlot$resolution <- newRes
            warning("Attempting to plot too many Hi-C pixels. Adjusting to ",
                    "a resolution of ", newRes, " BP.", call. = FALSE)
        } else {
            warning(hicPlot$resolution, " BP resolution detected in ",
            "dataframe. Attempting to plot too many Hi-C pixels. Please ",
            "read in data at a lower resolution before attempting to plot.", 
            call. = FALSE)
            hicPlot$resolution <- NA
        }
    
    }
    
    return(hicPlot)
}

## Define a function that reads in hic data for plotHic functions
# @param hic hic data argument
# @param hicPlot hic plot object
# @param norm normalization factor
# @param assembly genome assembly
# @param type matrix type
# @param quiet message quiet parameter
read_data <- function(hic, hicPlot, norm, assembly, type, quiet) {

    ## if .hic file, read in
    if (!("data.frame" %in% class(hic))) {
        if (!is.null(hicPlot$chromstart) & !is.null(hicPlot$chromend) &
            !is.na(hicPlot$resolution)) {
            readchromstart <- hicPlot$chromstart - hicPlot$resolution
            readchromend <- hicPlot$chromend + hicPlot$resolution
            readaltchromstart <- hicPlot$altchromstart - hicPlot$resolution
            readaltchromend <- hicPlot$altchromend + hicPlot$resolution
            hic <- suppressWarnings(readHic(
                file = hic, chrom = hicPlot$chrom,
                chromstart = readchromstart,
                chromend = readchromend,
                assembly = assembly,
                resolution = hicPlot$resolution,
                zrange = hicPlot$zrange,
                norm = norm,
                altchrom = hicPlot$altchrom,
                altchromstart = readaltchromstart,
                altchromend = readaltchromend,
                matrix = type
            ))
        } else {
            hic <- data.frame(matrix(nrow = 0, ncol = 3))
        }
    } else {
    
        if (!is.null(hicPlot$chromstart) & !is.null(hicPlot$chromend) &
            !is.na(hicPlot$resolution)) {
            if (!quiet) {
                message(
                    "Read in dataframe.  Assuming \'chrom\' in column1 ",
                    "and \'altchrom\' in column2. ",
                    hicPlot$resolution, " BP resolution detected."
                )
            }
            ## bbPlotHicRectangle specific warning for missing data
            if (hicPlot$chromstart < min(hic[, 1]) |
                hicPlot$chromend > max(hic[, 1])) {
                warning("`plotHicRectangle` requires additional data to",
                " plot a rectangular plot. Data is missing from input ",
                "dataframe for region and plot will be a trapezoid. To avoid ",
                "this missing data, call `plotHicRectangle` with full .hic",
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
# @param hicPlot hic plot object
set_zrange <- function(hic, hicPlot) {

    ## no zrange, only one value
    if (is.null(hicPlot$zrange) & length(unique(hic$counts)) == 1) {
        zrange <- c(unique(hic$counts), unique(hic$counts))
        hicPlot$zrange <- zrange
    }

    ## no zrange, multiple values
    if (is.null(hicPlot$zrange) & length(unique(hic$counts)) > 1) {
        if (grepl("log", hicPlot$colorTrans) == TRUE) {
            zrange <- c(0.0001, max(hic$counts))
        } else {
            zrange <- c(0, max(hic$counts))
        }

        hicPlot$zrange <- zrange
    }

    return(hicPlot)
}

## Define a function that parses an inherited half of a Hi-C plot
# @param hic hic plot object
inherit_half <- function(hic) {
    if (is(hic, "hicSquare")) {
        half <- hic$half
    } else if (is(hic, "hicTriangle") |
        is(hic, "hicRectangle")) {
        half <- "top"
    }

    return(half)
}
