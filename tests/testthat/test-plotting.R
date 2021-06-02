test_that("Plotting fxns produce expected plots", {
    library(BentoBoxData)
    save_png <- function(code, width = 400, height = 400) {
        path <- tempfile(fileext = ".png")
        png(path, width = width, height = height)
        on.exit(dev.off())
        code

        path
    }

    ## Hi-C
    data("bb_imrHicData")
    expect_snapshot_file(
        save_png(suppressMessages(bb_plotHicSquare(
            data = bb_imrHicData,
            chrom = "chr21",
            chromstart = 28000000,
            chromend = 30300000,
            resolution = 10000
        ))),
        "hic.png"
    )

    ## Signal
    data("bb_imrH3K27acData")
    expect_snapshot_file(
        save_png(suppressMessages(bb_plotSignal(
            data = bb_imrH3K27acData,
            chrom = "chr21",
            chromstart = 28000000,
            chromend = 30300000,
            range = c(0, 45)
        ))),
        "signal.png"
    )
})
