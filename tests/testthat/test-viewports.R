test_that("Viewport conversions", {
    pageCreate(width = 2, height = 2, default.units = "inches")
    testVP <- viewport(
        x = 1, y = 1,
        width = 1, height = 1,
        just = "center",
        default.units = "inches"
    )
    expect_equal(
        plotgardener:::vp_topLeft(testVP),
        list(unit(0.5, "inches"), unit(1.5, "inches"))
    )
    expect_equal(
        plotgardener:::vp_bottomLeft(testVP),
        list(unit(0.5, "inches"), unit(0.5, "inches"))
    )
    expect_equal(
        plotgardener:::vp_topRight(testVP),
        list(unit(1.5, "inches"), unit(1.5, "inches"))
    )
    expect_equal(
        plotgardener:::vp_bottomRight(testVP),
        list(unit(1.5, "inches"), unit(0.5, "inches"))
    )

    testVP <- viewport(
        x = 0.25, y = 0.25,
        width = 1.5, height = 1,
        just = c("left", "bottom"),
        default.units = "inches"
    )
    expect_equal(
        plotgardener:::adjust_vpCoords(testVP),
        list(unit(1, "inches"), unit(0.75, "inches"))
    )
})

test_that("Viewport order, naming, and numbering", {
    library(plotgardenerData)
    data("IMR90_DNAloops_pairs")
    pageCreate(width = 3, height = 5, default.units = "inches")
    arches <- suppressMessages(plotPairsArches(
        data = IMR90_DNAloops_pairs,
        chrom = "chr21", chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 0.5, y = 2.5, width = 2, height = 0.25,
        just = c("left", "top"), default.units = "inches",
        fill = "black", linecolor = "black", flip = TRUE
    ))
    expect_equal(plotgardener:::viewport_name(arches$grobs$vp), "arches1")
    expect_setequal(
        unlist(plotgardener:::current_viewports()),
        "arches1"
    )
    expect_equal(current.viewport()$name, "page")
    
    data("IMR90_ChIP_H3K27ac_signal")
    suppressMessages(plotSignal(
        data = IMR90_ChIP_H3K27ac_signal,
        chrom = "chr21", chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 0.5, y = 2.75, width = 2, height = 0.5,
        just = c("left", "top"), default.units = "inches"
    ))
    expect_setequal(
        unlist(plotgardener:::current_viewports()),
        c("arches1", "signal1_h")
    )
    suppressMessages(plotPairsArches(
        data = IMR90_DNAloops_pairs,
        chrom = "chr21", chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 0.5, y = 2.5, width = 2, height = 0.25,
        just = c("left", "top"), default.units = "inches",
        fill = "black", linecolor = "black", flip = TRUE
    ))
    expect_setequal(
        unlist(plotgardener:::current_viewports()),
        c("arches1", "signal1_h", "arches2")
    )

    pagePlotRemove(plot = arches)
    expect_setequal(
        unlist(plotgardener:::current_viewports()),
        c("signal1_h", "arches2")
    )
    
})

test_that("Below-y coordinate calculation", {
    library(plotgardenerData)
    data("IMR90_HiC_10kb")
    pageCreate(width = 3, height = 5, default.units = "inches")
    suppressMessages(plotHicSquare(
        data = IMR90_HiC_10kb,
        chrom = "chr21",
        chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top"), default.units = "inches"
    ))
    expect_equal(plotgardener:::plot_belowY("0b"), unit(2.5, "inches"))
})

test_that("draw parameter and pagePlotPlace", {
    library(plotgardenerData)
    data("IMR90_ChIP_H3K27ac_signal")
    pageCreate(width = 3, height = 3, default.units = "inches")

    expect_error(signalPlot <- plotSignal(
        data = IMR90_ChIP_H3K27ac_signal,
        chrom = "chr21", chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 0.25,
        draw = TRUE
    ))

    signalPlot <- suppressMessages(plotSignal(
        data = IMR90_ChIP_H3K27ac_signal,
        chrom = "chr21", chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 0.25,
        draw = FALSE
    ))
    expect_equal(
        unlist(plotgardener:::current_viewports()),
        NULL
    )

    signalPlot <- suppressMessages(pagePlotPlace(
        plot = signalPlot,
        x = 0.5, y = 0.5, width = 2, height = 1,
        just = c("left", "top"), default.units = "inches"
    ))
    expect_setequal(
        unlist(plotgardener:::current_viewports()),
        c("signal1_h")
    )
    expect_equal(signalPlot$x, unit(0.5, "inches"))
})

test_that("page unit conversions", {
    
    pageCreate(width = 3, height = 3, default.units = "inches")
    # Same units but convert to proper y-coordinate
    testObject <- list("x" = unit(0.5, "inches"),
                   "y" = unit(0.5, "inches"),
                   "width" = unit(1, "inches"),
                   "height" = unit(1, "inches"))
    expect_equal(plotgardener:::convert_page(testObject),
                 list("x" = unit(0.5, "inches"),
                      "y" = unit(2.5, "inches"),
                      "width" = unit(1, "inches"),
                      "height" = unit(1, "inches")))
    
    # Different units and proper y-coordinate conversion
    testObject2 <- list("x" = unit(0.5, "npc"),
                        "y" = unit(0.5, "npc"),
                        "width" =  unit(1, "npc"),
                        "height" = unit(1, "npc"))
    expect_equal(plotgardener:::convert_page(testObject2),
                 list("x" = unit(1.5, "inches"),
                      "y" = unit(1.5, "inches"),
                      "width" = unit(3, "inches"),
                      "height" = unit(3, "inches")))
})

test_that("annotation viewports", {
    
    ## Not adding annotation viewports with addViewport
    library(plotgardenerData)
    data("IMR90_HiC_10kb")
    pageCreate(width = 3, height = 5, default.units = "inches")
    suppressMessages(hicPlot <- plotHicSquare(
        data = IMR90_HiC_10kb,
        chrom = "chr21",
        chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 0.5, y = 0.5, width = 2, height = 2,
        just = c("left", "top"), default.units = "inches"
    ))
    suppressMessages(annoHeatmapLegend(plot = hicPlot, x = 2.6, y = 0.5, 
                                            height = 0.5, width = 0.1))
    expect_equal(unlist(get("pg_vpTree", envir = pgEnv)),
                 c("hicSquare1"))
    
    ## Getting appropriate annotation viewports for genes, hic triangle/rects
    
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    library("org.Hs.eg.db")

    ## Set genomic coordinates
    paramssmall <- pgParams(
        chrom = "chr8",
        chromstart = 1, chromend = 3000000,
        assembly = "hg19", width = 7
    )
    
    pageCreate(width = 7.5, height = 3.5, default.units = "inches")
    genesPlot <- suppressMessages(plotGenes(
            params = paramssmall,
            geneHighlights = data.frame(
                "gene" = c("DLGAP2"),
                "color" = c("#225EA8")
            ),
            geneBackground = "grey",
            x = 0.25, y = 2.25, height = 0.75,
            just = c("left", "top"), default.units = "inches"
        ))
    expect_equal(plotgardener:::getAnnoViewport(genesPlot),
                 genesPlot$grobs$children$background$vp)
    

    pageCreate(width = 4, height = 2.5, default.units = "inches")
    hicPlot <- suppressMessages(plotHicTriangle(
        data = IMR90_HiC_10kb, resolution = 10000,
        zrange = c(0, 70),
        chrom = "chr21",
        chromstart = 28000000, chromend = 30300000,
        assembly = "hg19",
        x = 2, y = 0.5, width = 3, height = 1.5,
        just = "top", default.units = "inches"
    ))
    expect_equal(plotgardener:::getAnnoViewport(hicPlot),
                 hicPlot$outsideVP)
})
