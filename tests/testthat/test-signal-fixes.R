test_that("plotSignal baseline is drawn last and can be disabled", {
    signalData <- data.frame(
        chrom = "chr1",
        start = seq(1, 901, by = 100),
        end = seq(100, 1000, by = 100),
        score = rep(10, 10)
    )

    pageCreate(width = 4, height = 3, default.units = "inches",
               showGuides = FALSE)
    withBaseline <- plotSignal(
        signalData, chrom = "chr1", chromstart = 1, chromend = 1000,
        assembly = "hg38", fill = "steelblue", binSize = 100,
        x = 0.5, y = 0.5, width = 3, height = 1,
        draw = FALSE, baseline = TRUE
    )
    withoutBaseline <- plotSignal(
        signalData, chrom = "chr1", chromstart = 1, chromend = 1000,
        assembly = "hg38", fill = "steelblue", binSize = 100,
        x = 0.5, y = 1.6, width = 3, height = 1,
        draw = FALSE, baseline = FALSE
    )

    expect_equal(
        length(withBaseline$grobs$children),
        length(withoutBaseline$grobs$children) + 1
    )
})

test_that("plotSignal label position is configurable", {
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
        x = 0.5, y = 0.5, width = 3, height = 1,
        draw = FALSE, label = "custom", label.position = "bottom-right"
    )
    expect_true(any(vapply(
        sig$grobs$children,
        function(g) identical(g$label, "custom"),
        logical(1)
    )))
})

test_that("plotMultiSignal honors a user-supplied range", {
    signalData <- data.frame(
        chrom = "chr1",
        start = seq(1, 901, by = 100),
        end = seq(100, 1000, by = 100),
        score = rep(10, 10)
    )
    pageCreate(width = 4, height = 3, default.units = "inches",
               showGuides = FALSE)
    sigs <- plotMultiSignal(
        list(signalData, signalData),
        chrom = "chr1", chromstart = 1, chromend = 1000,
        assembly = "hg38", range = c(0, 250), binSize = 100,
        x = 0.5, y = 0.5, width = 3, height = 2,
        draw = FALSE, label = c("A", "B")
    )
    expect_equal(sigs[[1]]$range, c(0, 250))
    expect_equal(sigs[[2]]$range, c(0, 250))
})
