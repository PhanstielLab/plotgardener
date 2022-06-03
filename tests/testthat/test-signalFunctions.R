test_that("signalRanges", {
    library(plotgardenerData)
    data("IMR90_ChIP_H3K27ac_signal")
    
    region <- pgParams(
        chrom = "chr21",
        chromstart = 28000000, chromend = 30300000,
        assembly = "hg19")
    
    sig1 <- suppressMessages(plotSignal(
        data = IMR90_ChIP_H3K27ac_signal, params = region, range = c(0, 45)))
    
    expect_equal(sig1$range, c(0, 45))
    
})

test_that("plotSignalMessages", {
    
    library(plotgardenerData)
    data("IMR90_ChIP_H3K27ac_signal")
    
    region <- pgParams(
        chrom = "chr21",
        chromstart = 28000000, chromend = 30300000,
        assembly = "hg19")
    
    # Positive data should have no warnings/errors
    expect_message(plotSignal(data = IMR90_ChIP_H3K27ac_signal,
                              params = region))
    expect_invisible(plotSignal(data = IMR90_ChIP_H3K27ac_signal,
                                params = region))
    
})