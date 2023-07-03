library(withr)
library(pdftools)
library(plotgardenerData)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(AnnotationHub)
library(org.Hs.eg.db)


test_that("plotGenes_pdf", {
    temp_pdf <- file.path(tempdir(), "testGenes.pdf")
    on.exit(unlink(temp_pdf))

    graphics.off()
    pdf(temp_pdf, width = 7, height = 3.5)
    suppressMessages(plotGenes(chrom = "chr8",
             chromstart = 1,
              chromend = 3000000,
              assembly = "hg19"))
    dev.off()
    
    expect_equal(pdf_length(temp_pdf), 1)
    
})

test_that("plotHicRectangle_pdf", {

    data("IMR90_HiC_10kb")
    temp_pdf <- tempfile("testHicRectangle", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 6, height = 3.5)
    suppressMessages(plotHicRectangle(data = IMR90_HiC_10kb, resolution = 10000,
                     zrange = c(0, 70),
                     chrom = "chr21",
                     chromstart = 28950000,
                     chromend = 29800000,
                     assembly = "hg19"))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})

test_that("plotHicSquare_pdf", {

    data("IMR90_HiC_10kb")
    temp_pdf <- tempfile("testHicSquare", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 3, height = 3)
    suppressMessages(plotHicSquare(data = IMR90_HiC_10kb, resolution = 10000,
                  zrange = c(0, 70),
                  chrom = "chr21",
                  chromstart = 28000000,
                  chromend = 30300000,
                  assembly = "hg19"))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})

test_that("plotHicTriangle_pdf", {

    data("IMR90_HiC_10kb")
    temp_pdf <- tempfile("testHicTriangle", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 4, height = 2.5)
    suppressMessages(plotHicTriangle(data = IMR90_HiC_10kb, resolution = 10000,
                    zrange = c(0, 70),
                    chrom = "chr21",
                    chromstart = 28000000,
                    chromend = 30300000,
                    assembly = "hg19"))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})


test_that("plotIdeogram_pdf", {

    temp_pdf <- tempfile("testIdeogram", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 4.5, height = 1)
    suppressMessages(plotIdeogram(chrom = "chr2", assembly = "hg19"))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})

test_that("plotLegend_pdf", {

    temp_pdf <- tempfile("testLegend", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 2, height = 2)
    suppressMessages(plotLegend(legend = c("- strand", "+ strand"),
               fill = c("steel blue", "light salmon"),
               border = FALSE))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)
})

test_that("plotManhattan_pdf", {

    data("hg19_insulin_GWAS")

    temp_pdf <- tempfile("testManhattan", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 7.5, height = 4.5)
    suppressMessages(plotManhattan(data = hg19_insulin_GWAS, assembly = "hg19",
                  fill = c("grey", "#37a7db"),
                  sigLine = TRUE, col = "grey", range = c(0, 14)))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})

test_that("plotPairs_pdf", {

    data("IMR90_DNAloops_pairs")

    temp_pdf <- tempfile("testPairs", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 7.5, height = 2)
    suppressWarnings(plotPairs(data = IMR90_DNAloops_pairs,
              chrom = "chr21",
              chromstart = 27900000, chromend = 30700000,
              assembly = "hg19"))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})

test_that("plotPairsArches_pdf", {

    data("IMR90_DNAloops_pairs")

    temp_pdf <- tempfile("testPairsArches", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 7.5, height = 2)
    suppressWarnings(plotPairsArches(data = IMR90_DNAloops_pairs,
                    chrom = "chr21",
                    chromstart = 27900000, chromend = 30700000,
                    assembly = "hg19"))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})

test_that("plotRanges_pdf", {

    data("IMR90_ChIP_CTCF_reads")

    temp_pdf <- tempfile("testRanges", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 7.5, height = 5)
    suppressWarnings(plotRanges(data = IMR90_ChIP_CTCF_reads,
               chrom = "chr21",
               chromstart = 29073000, chromend = 29074000,
               assembly = "hg19"))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})

test_that("plotSignal_pdf", {

    data("IMR90_ChIP_H3K27ac_signal")

    temp_pdf <- tempfile("testSignal", fileext = ".pdf")
    on.exit(unlink(temp_pdf))
    graphics.off()
    pdf(temp_pdf, width = 7.5, height = 2)
    suppressMessages(plotSignal(data = IMR90_ChIP_H3K27ac_signal,
               chrom = "chr21",
               chromstart = 28000000, chromend = 30300000,
               assembly = "hg19",
               range = c(0, 45)))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})

test_that("plotTranscripts_pdf", {

    temp_pdf <- tempfile("testTranscripts", fileext = ".pdf")
    graphics.off()
    pdf(temp_pdf, width = 7.5, height = 5)
    on.exit(unlink(temp_pdf))
    suppressWarnings(plotTranscripts(chrom = "chr8",
                    chromstart = 1000000,
                    chromend = 2000000,
                    assembly = "hg19"))
    dev.off()

    expect_equal(pdf_length(temp_pdf), 1)

})