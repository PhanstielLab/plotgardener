test_that("plotHicSquare rasterizes when requested", {
    pairData <- expand.grid(x = seq(0, 90, 10), y = seq(0, 90, 10))
    pairData <- pairData[pairData$x <= pairData$y, ]
    pairData$counts <- 1

    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = FALSE)
    hic <- plotHicSquare(
        data = pairData, resolution = 10, zrange = c(0, 2),
        chrom = "chr1", chromstart = 0, chromend = 100,
        assembly = "hg38", raster = "always",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top"), draw = FALSE
    )
    expect_s3_class(hic$grobs$children[[2]], "rastergrob")

    hicVector <- plotHicSquare(
        data = pairData, resolution = 10, zrange = c(0, 2),
        chrom = "chr1", chromstart = 0, chromend = 100,
        assembly = "hg38", raster = "never",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top"), draw = FALSE
    )
    expect_s3_class(hicVector$grobs$children[[2]], "rect")
})

test_that("plotHicRectangle rasterizes when requested", {
    pairData <- expand.grid(x = seq(0, 90, 10), y = seq(0, 90, 10))
    pairData <- pairData[pairData$x <= pairData$y, ]
    pairData$counts <- 1

    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = FALSE)
    hic <- suppressWarnings(plotHicRectangle(
        data = pairData, resolution = 10, zrange = c(0, 2),
        chrom = "chr1", chromstart = 0, chromend = 100,
        assembly = "hg38", raster = "always",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top"), draw = FALSE
    ))
    expect_s3_class(hic$grobs$children[[2]], "rastergrob")
})

test_that("plotHicTriangle rasterizes regular pixels", {
    pairData <- expand.grid(x = seq(0, 90, 10), y = seq(0, 90, 10))
    pairData <- pairData[pairData$x <= pairData$y, ]
    pairData$counts <- 1

    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = FALSE)
    hic <- plotHicTriangle(
        data = pairData, resolution = 10, zrange = c(0, 2),
        chrom = "chr1", chromstart = 0, chromend = 100,
        assembly = "hg38", raster = "always",
        x = 0.5, y = 0.5, width = 2, height = 1,
        just = c("left", "top"), draw = FALSE
    )
    expect_true(any(vapply(
        hic$grobs$children,
        function(g) inherits(g, "rastergrob"),
        logical(1)
    )))
})
