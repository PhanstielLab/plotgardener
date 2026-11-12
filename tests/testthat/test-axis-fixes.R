test_that("annoYaxis supports tick and axis-line styling", {
    signalData <- data.frame(
        chrom = "chr1",
        start = seq(1, 901, by = 100),
        end = seq(100, 1000, by = 100),
        score = rep(10, 10)
    )
    pageCreate(width = 4, height = 3, default.units = "inches",
               showGuides = FALSE)
    sig <- plotSignal(
        signalData, chrom = "chr1", chromstart = 1, chromend = 1000,
        assembly = "hg38", binSize = 100,
        x = 0.5, y = 0.5, width = 3, height = 1, draw = FALSE
    )
    yax <- annoYaxis(
        plot = sig, at = c(0, 10), axisLine = TRUE,
        tcl = 1, label.offset = 0.25,
        axisLine.lwd = 2, axisLine.col = "blue",
        tick.lwd = 3, tick.col = "red"
    )

    yGrob <- yax$grobs$children[[1]]
    expect_equal(yGrob$children$ticks$x1, unit(-1, "lines"))
    expect_equal(yGrob$children$labels$x, unit(-1.25, "lines"))
    expect_equal(yGrob$children$major$gp$lwd, 2)
    expect_equal(yGrob$children$major$gp$col, "blue")
    expect_equal(yGrob$children$ticks$gp$lwd, 3)
    expect_equal(yGrob$children$ticks$gp$col, "red")
})

test_that("annoXaxis axisLine = FALSE removes the axis line without grid.remove", {
    signalData <- data.frame(
        chrom = "chr1",
        start = seq(1, 901, by = 100),
        end = seq(100, 1000, by = 100),
        score = rep(10, 10)
    )
    pageCreate(width = 4, height = 3, default.units = "inches",
               showGuides = FALSE)
    sig <- plotSignal(
        signalData, chrom = "chr1", chromstart = 1, chromend = 1000,
        assembly = "hg38", binSize = 100,
        x = 0.5, y = 0.5, width = 3, height = 1, draw = FALSE
    )
    xax <- annoXaxis(plot = sig, at = c(1, 1000), axisLine = FALSE)
    xGrob <- xax$grobs$children[[1]]
    expect_false("major" %in% names(xGrob$children))
})

test_that("annoYaxis defaults use shorter ticks and thinner lines", {
    signalData <- data.frame(
        chrom = "chr1",
        start = seq(1, 901, by = 100),
        end = seq(100, 1000, by = 100),
        score = rep(10, 10)
    )
    pageCreate(width = 4, height = 3, default.units = "inches",
               showGuides = FALSE)
    sig <- plotSignal(
        signalData, chrom = "chr1", chromstart = 1, chromend = 1000,
        assembly = "hg38", binSize = 100,
        x = 0.5, y = 0.5, width = 3, height = 1, draw = FALSE
    )
    yax <- annoYaxis(plot = sig, at = c(0, 10), axisLine = TRUE)
    yGrob <- yax$grobs$children[[1]]
    expect_equal(yGrob$children$ticks$x1, unit(-0.25, "lines"))
    expect_equal(yGrob$children$labels$x, unit(-0.6, "lines"))
    expect_equal(yGrob$children$major$gp$lwd, 0.5)
    expect_equal(yGrob$children$ticks$gp$lwd, 0.5)
})
