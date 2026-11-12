test_that("plotHicSquare does not mirror interchromosomal matrices", {
    pairData <- data.frame(
        x = rep(seq(100, 300, by = 100), each = 3),
        y = rep(seq(500, 700, by = 100), 3),
        counts = rep(1, 9)
    )

    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = FALSE)
    hic <- plotHicSquare(
        data = pairData, resolution = 100, zrange = c(0, 2),
        chrom = "chr4", chromstart = 0, chromend = 400,
        altchrom = "chr14", altchromstart = 400, altchromend = 800,
        assembly = "hg38", half = "top",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top"), draw = FALSE
    )

    expect_equal(hic$xchrom, "chr4")
    expect_equal(hic$ychrom, "chr14")
    expect_equal(hic$xchromstart, 0)
    expect_equal(hic$ychromstart, 400)

    ## No symmetrization: one rectangle per input pixel
    rectGrob <- hic$grobs$children[[2]]
    expect_equal(length(rectGrob$x), nrow(pairData))
})

test_that("annoGenomeLabel uses explicit interchromosomal axis metadata", {
    pairData <- data.frame(
        x = rep(seq(100, 300, by = 100), each = 3),
        y = rep(seq(500, 700, by = 100), 3),
        counts = rep(1, 9)
    )

    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = FALSE)
    hic <- plotHicSquare(
        data = pairData, resolution = 100, zrange = c(0, 2),
        chrom = "chr4", chromstart = 0, chromend = 400,
        altchrom = "chr14", altchromstart = 400, altchromend = 800,
        assembly = "hg38", half = "top",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top"), draw = FALSE
    )

    labels <- annoGenomeLabel(
        hic, axis = "y", scale = "Kb",
        x = 0.45, y = 0.5, just = c("right", "top")
    )
    expect_s3_class(labels, "genomeLabel")
})

test_that("plotHicSquare fixed aspect uses equal basepairs per axis", {
    pairData <- data.frame(
        x = rep(seq(100, 300, by = 100), each = 3),
        y = rep(seq(500, 900, by = 200), 3),
        counts = rep(1, 9)
    )

    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = FALSE)
    hic <- plotHicSquare(
        data = pairData, resolution = 100, zrange = c(0, 2),
        chrom = "chr4", chromstart = 0, chromend = 400,
        altchrom = "chr14", altchromstart = 400, altchromend = 1000,
        assembly = "hg38", half = "top", aspect = "fixed",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top"), draw = FALSE
    )

    expect_equal(diff(range(hic$xaxisRange)),
                 diff(range(hic$yaxisRange)))
    expect_equal(diff(range(hic$xaxisRange)), 600)

    xlab <- annoGenomeLabel(hic, axis = "x", x = 0.5, y = 2.8,
                            just = c("left", "top"))
    ylab <- annoGenomeLabel(hic, axis = "y", x = 0.45, y = 0.5,
                            just = c("right", "top"))
    expect_s3_class(xlab, "genomeLabel")
    expect_s3_class(ylab, "genomeLabel")
})
