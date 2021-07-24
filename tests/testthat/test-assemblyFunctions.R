test_that("chromDataAgreement", {
    library(BentoBoxData)
    ## Ranges
    data("IMR90_ChIP_CTCF_reads")
    ## match
    expect_silent(BentoBox:::chromDataAgreement(data = IMR90_ChIP_CTCF_reads,
                                                chrom = "chr21",
                                                type = "ranges"))
    ## don't match
    expect_warning(BentoBox:::chromDataAgreement(data = IMR90_ChIP_CTCF_reads,
                                                chrom = "21",
                                                type = "ranges"))
    ## Pairs
    data("IMR90_DNAloops_pairs")
    ## match
    expect_silent(BentoBox:::chromDataAgreement(data = IMR90_DNAloops_pairs,
                                                chrom = "chr21",
                                                type = "pairs"))
    ## don't match
    expect_warning(BentoBox:::chromDataAgreement(data = IMR90_DNAloops_pairs,
                                                chrom = "21",
                                                type = "pairs"))
})

test_that("genomicScale", {
    ## Doesn't have TxDb and has c(0, 1) xscale
    testObject <- list("chromstart" = NULL,
                       "chromend" = NULL,
                       "assembly" = BentoBox:::parse_bbAssembly("mm9"))
    testObjectInternal <- list()
    expect_warning(BentoBox:::genomicScale(object = testObject,
                                           objectInternal = testObjectInternal,
                                           plotType = "genes"))
    scale <- suppressWarnings(BentoBox:::genomicScale(object = testObject,
                                     objectInternal = testObjectInternal,
                                     plotType = "genes"))[[2]]$xscale
    expect_equal(scale, c(0, 1))
    
    ## Has TxDb and full chromosome xscale
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    testObject <- list("chrom" = "chr1",
                    "chromstart" = NULL,
                    "chromend" = NULL,
                    "assembly" = BentoBox:::parse_bbAssembly("hg19"))
    testObjectInternal <- list()
    scale <- BentoBox:::genomicScale(object = testObject,
                                     objectInternal = testObjectInternal,
                                     plotType = "genes")[[2]]$xscale
    expect_equal(scale, c(1, 249250621))
    
    ## Throws wrong chromosome warning
    testObject <- list("chrom" = "chrHello",
                       "chromstart" = NULL,
                       "chromend" = NULL,
                       "assembly" = BentoBox:::parse_bbAssembly("hg19"))
    expect_warning(BentoBox:::genomicScale(object = testObject,
                                           objectInternal = testObjectInternal,
                                           plotType = "genes"))
    scale2 <- suppressWarnings(BentoBox:::genomicScale(object = testObject,
                                        objectInternal = testObjectInternal,
                                        plotType = "genes"))[[2]]$xscale
    expect_equal(scale2, c(0, 1))
})

test_that("geneData", {
    
    ## TxDb and OrgDb gets data
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    library("org.Hs.eg.db")
    testObject <- list("chrom" = "chr1",
                        "assembly" = BentoBox:::parse_bbAssembly("hg19"))
    testObjectInternal <- list()
    
    expect_silent(BentoBox:::geneData(object = testObject,
                                      objectInternal = testObjectInternal))
})