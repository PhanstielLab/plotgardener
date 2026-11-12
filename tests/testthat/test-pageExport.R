test_that("pageExport writes a single-page PDF without dev.off", {
    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = TRUE)
    plotRect(x = 0.5, y = 0.5, width = 2, height = 2,
             fill = "green", default.units = "inches")

    out <- withr::local_tempfile(fileext = ".pdf")
    pageExport(out)
    expect_true(file.exists(out))
    expect_gt(file.info(out)$size, 0)

    raw <- readBin(out, what = "raw", n = file.info(out)$size)
    pages <- length(grepRaw(charToRaw("/Type /Page /Parent"), raw,
                            all = TRUE))
    expect_equal(pages, 1)
})

test_that("pageExport supports raster output and guide hiding", {
    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = TRUE)
    plotRect(x = 0.5, y = 0.5, width = 2, height = 2,
             fill = "green", default.units = "inches")

    out <- withr::local_tempfile(fileext = ".png")
    pageExport(out, dpi = 100, guides = FALSE)
    expect_true(file.exists(out))
    expect_gt(file.info(out)$size, 0)
})

test_that("pageGuideHide and pagePlotRemove are deferred to export", {
    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = TRUE)
    r <- plotRect(x = 0.5, y = 0.5, width = 1, height = 1,
                  fill = "green", default.units = "inches")
    pageGuideHide(redraw = FALSE)
    pagePlotRemove(r, redraw = FALSE)

    expect_true(get("guides_hidden", envir = plotgardener:::pgEnv))
    expect_true(r$grobs$name %in%
                get("removed_plots", envir = plotgardener:::pgEnv))

    out <- withr::local_tempfile(fileext = ".pdf")
    pageExport(out)
    raw <- readBin(out, what = "raw", n = file.info(out)$size)
    pages <- length(grepRaw(charToRaw("/Type /Page /Parent"), raw,
                            all = TRUE))
    expect_equal(pages, 1)
})
