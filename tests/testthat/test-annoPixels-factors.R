test_that("annoPixels handles factor chromosome levels from different sources", {
    pairData <- expand.grid(x = seq(0, 90, 10), y = seq(0, 90, 10))
    pairData <- pairData[pairData$x <= pairData$y, ]
    pairData$counts <- 1

    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = FALSE)
    hic <- suppressWarnings(plotHicRectangle(
        data = pairData, resolution = 10, zrange = c(0, 2),
        chrom = "chr1", chromstart = 0, chromend = 100,
        assembly = "hg38",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top")
    ))

    loops <- data.frame(
        chrom1 = factor(c("chr1", "chr2"), levels = c("chr1", "chr2")),
        start1 = c(10, 20),
        end1 = c(20, 30),
        chrom2 = factor(c("chr1", "chr2"), levels = c("chr1", "chr2")),
        start2 = c(10, 20),
        end2 = c(20, 30)
    )
    expect_silent(suppressMessages(
        annoPixels(plot = hic, data = loops)
    ))
})
