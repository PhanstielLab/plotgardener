test_that("signalRanges", {
    library(plotgardenerData)
    data("IMR90_ChIP_H3K27ac_signal")
    data("GM12878_ChIP_H3K27ac_signal")
    GM12878_ChIP_H3K27ac_signal[,"score"] <- 
        GM12878_ChIP_H3K27ac_signal[,"score"] * -1
    
    region <- pgParams(
            chrom = "chr21",
            chromstart = 28000000, chromend = 30300000,
            assembly = "hg19")
    
    sig1 <- suppressMessages(plotSignal(
        data = IMR90_ChIP_H3K27ac_signal, params = region, range = c(0, 45)))

    expect_equal(sig1$range, c(0, 45))

    sig2 <- suppressMessages(plotSignal(
        data = list(IMR90_ChIP_H3K27ac_signal, 
                    GM12878_ChIP_H3K27ac_signal), params = region, 
        range = c(-45, 45)
    ))
    expect_equal(sig2$range, c(-45, 45))

})

test_that("plotSignalMessages", {
    
    library(plotgardenerData)
    data("IMR90_ChIP_H3K27ac_signal")
    data("GM12878_ChIP_H3K27ac_signal")
    GM12878_ChIP_H3K27ac_signal[,"score"] <- 
                                    GM12878_ChIP_H3K27ac_signal[,"score"] * -1
    
    region <- pgParams(
        chrom = "chr21",
        chromstart = 28000000, chromend = 30300000,
        assembly = "hg19")
    
    # Positive data should have no warnings/errors
    expect_message(plotSignal(data = IMR90_ChIP_H3K27ac_signal,
                              params = region))
    expect_invisible(plotSignal(data = IMR90_ChIP_H3K27ac_signal,
                              params = region))
    
    # Negative data with negData = TRUE should have no warnings or errors
    expect_message(plotSignal(data = GM12878_ChIP_H3K27ac_signal,
                              params = region, negData = TRUE))
    
    # Negative data with negData = FALSE should throw warning
    expect_warning(plotSignal(data = GM12878_ChIP_H3K27ac_signal,
                              params = region))
    
    # Positive and negative data in proper order should have no warnings/errors
    expect_message(plotSignal(data = list(IMR90_ChIP_H3K27ac_signal,
                                          GM12878_ChIP_H3K27ac_signal),
                              params = region))
    
    # Positive and negative data in proper order should error
    expect_error(plotSignal(data = list(GM12878_ChIP_H3K27ac_signal,
                                        IMR90_ChIP_H3K27ac_signal),
                              params = region))
    
})