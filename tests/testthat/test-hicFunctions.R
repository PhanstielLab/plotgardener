test_that("adjust_resolution", {
    library(plotgardenerData)
    data("IMR90_HiC_10kb")
    ## Detect appropriate resolution from data
    expect_equal(plotgardener:::adjust_resolution(hic = IMR90_HiC_10kb,
                                              hicPlot = NULL)$resolution,
                 10000)
    ## Grab a resolution in a file
    hicFile <- system.file("extdata/test_chr22.hic", package="plotgardenerData")
    hic_plot <- list("chromstart" = 20000000,
                     "chromend" = 47500000)
    expect_equal(plotgardener:::adjust_resolution(hic = hicFile,
                                              hicPlot = hic_plot)$resolution,
                 100000)
    
    ## Grab a close resolution to one in a file
    hic_plot <- list("chromstart" = 20000000,
                     "chromend" = 40000000)
    expect_equal(plotgardener:::adjust_resolution(hic = hicFile,
                                              hicPlot = hic_plot)$resolution,
                 100000)
    
})

test_that("hic_limit", {
    library(plotgardenerData)
    hicFile <- system.file("extdata/test_chr22.hic", package="plotgardenerData")
    
    ## Appropriate number of bins
    hic_plot <- list("chromstart" = 20000000,
                     "chromend" = 47500000,
                     "resolution" = 100000)
    expect_silent(plotgardener:::hic_limit(hic = hicFile,
                                       hicPlot = hic_plot))
    
    ## Hic file and readjusts bin number
    hic_plot <- list("chromstart" = 20000000,
                     "chromend" = 47500000,
                     "resolution" = 5000)
    expect_warning(plotgardener:::hic_limit(hic = hicFile,
                                        hicPlot = hic_plot))
    expect_equal(suppressWarnings(plotgardener:::hic_limit(hic = hicFile,
                                      hicPlot = hic_plot)$resolution),
                 100000)
    
    ## Data frame NA resolution
    data("IMR90_HiC_10kb")
    expect_warning(plotgardener:::hic_limit(hic = IMR90_HiC_10kb,
                                        hicPlot = hic_plot))
    expect_equal(suppressWarnings(plotgardener:::hic_limit(hic = IMR90_HiC_10kb,
                                      hicPlot = hic_plot)$resolution),
                 NA)
})

test_that("set_zrange", {
    library(plotgardenerData)
    data("IMR90_HiC_10kb")
    ## already has a zrange
    hic_plot <- list("zrange" = c(0, 10))
    
    expect_equal(plotgardener:::set_zrange(IMR90_HiC_10kb, 
                                       hicPlot = hic_plot)$zrange,
                 c(0, 10))
    
    
    ## no zrange, multiple values
    hic_plot <- list("zrange" = NULL, "colorTrans" = "linear")
    expect_equal(plotgardener:::set_zrange(IMR90_HiC_10kb, 
                                       hicPlot = hic_plot)$zrange,
                 c(0, 70))
    
    ## no zrange, only one value
    IMR90_HiC_10kb$counts <- 1
    expect_equal(plotgardener:::set_zrange(IMR90_HiC_10kb, 
                                       hicPlot = hic_plot)$zrange,
                 c(1, 1))
})

test_that("inherit_half", {
    
    ## square hic half , top
    library(plotgardenerData)
    data("IMR90_HiC_10kb")
    pageCreate(width = 3, height = 3, default.units = "inches")
    hicPlot <- suppressMessages(plotHicSquare(
            data = IMR90_HiC_10kb, resolution = 10000,
            zrange = c(0, 70),
            chrom = "chr21",
            chromstart = 28000000, chromend = 30300000,
            assembly = "hg19",
            half = "top",
            x = 0.5, y = 0.5, width = 2, height = 2,
            just = c("left", "top"),
            default.units = "inches"
        ))
    expect_equal(plotgardener:::inherit_half(hicPlot),
                 "top")
    
    ## triangle hic half, always top
    pageCreate(width = 4, height = 2.5, default.units = "inches")

    ## Plot and place triangle Hi-C plot
    hicPlot <- suppressMessages(plotHicTriangle(
        data = IMR90_HiC_10kb, resolution = 10000,
        zrange = c(0, 70),
        chrom = "chr21",
        chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 2, y = 0.5, width = 3, height = 1.5,
        just = "top", default.units = "inches"
    ))
    hicPlot$half <- "testing"
    expect_equal(plotgardener:::inherit_half(hicPlot),
                 "top")
})
