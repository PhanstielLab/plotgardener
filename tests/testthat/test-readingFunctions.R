test_that("read_rangeData", {
    library(GenomicRanges)
    library(plotgardenerData)
    ## GRanges reading
    gr <- GRanges(seqnames = "chr1", strand = c("+", "-", "+"),
                  ranges = IRanges(start = c(1,3,5), width = 3))
    values(gr) <- DataFrame(score = c(0.1, 0.5, 0.3))
    
    expectedgr <- as.data.frame(gr)
    colnames(expectedgr)[1:3] <- c("chrom", "start", "end")
    
    expect_equal(plotgardener:::read_rangeData(data = gr,
                              assembly = "hg19"),
                 expectedgr)
    
    ## Errors for invalid column types
    expectedgr$start <- as.character(expectedgr$start)
    expectedgr$end <- as.character(expectedgr$end)
    expect_error(plotgardener:::read_rangeData(data = expectedgr,
                              assembly = "hg19"))
    
})

test_that("read_pairedData", {
    ## GInteractions reading
    library(GenomicRanges)
    library(InteractionSet)
    all.regions <- GRanges("chrA", IRanges(0:9*10+1, 1:10*10))
    index.1 <- c(1,2,3)
    index.2 <- c(4,5,6)
    region.1 <- all.regions[index.1]
    region.2 <- all.regions[index.2]
    gi <- GInteractions(region.1, region.2)
    
    expectedgi <- as.data.frame(gi)[,c(1,2,3,6,7,8,4,5,9,10)]
    
    colnames(expectedgi)[1:6] <- c("chrom1", "start1", "end1",
                                   "chrom2", "start2", "end2")
    expect_equal(plotgardener:::read_pairedData(data = gi,
                                           assembly = "hg19"),
                 expectedgi)
    
    ## Warning for out of order anchors
    index.3 <- c(1,6,8)
    index.4 <- c(2,4,7)
    region.3 <- all.regions[index.3]
    region.4 <- all.regions[index.4]
    gi2 <- GInteractions(region.3, region.4)
    expect_warning(plotgardener:::read_pairedData(data = gi2,
                                                  assembly = "hg19"))
    
    ## Errors for invalid column types
    expectedgi$start1 <- as.character(expectedgi$start1)
    expectedgi$end2 <- as.character(expectedgi$end2)
    expect_error(plotgardener:::read_pairedData(data = expectedgi,
                                            assembly = "hg19"))
    
})

test_that("checkAssemblyMatch", {
    library("TxDb.Hsapiens.UCSC.hg19.knownGene")
    ## warning for invalid matching
    tx_db <- TxDb.Hsapiens.UCSC.hg19.knownGene
    expect_warning(plotgardener:::checkAssemblyMatch(data = tx_db,
            assembly = plotgardener:::parseAssembly("hg38")))
})

test_that("readHic", {
    library("plotgardenerData")
    
    hicFile <- system.file("extdata/test_chr22.hic", package="plotgardenerData")
    
    ## Check that Hi-C file reads in with correct dimension
    hicData <- readHic(file = hicFile, chrom = "22",
                       chromstart = 20000000, chromend = 20300000,
                       assembly = "hg19",
                       resolution = 100000) 
    expect_equal(nrow(hicData), 10)
    expect_equal(min(hicData[,1]), 20000000)
    expect_equal(max(hicData[,1]), 20300000)
})

test_that("readCool", {
    skip_if_offline()

    ## .cool file
    coolFile <- file.path(tempdir(), "Rao2014-IMR90-MboI-allreps-filtered.1000kb.cool")
    download.file(url = "https://usgs2.osn.mghpcc.org/cooler01/examples/hg19/Rao2014-IMR90-MboI-allreps-filtered.1000kb.cool",
                            destfile = coolFile, mode = "wb")
    on.exit(unlink(coolFile))
    
    # File type
    expect_equal(.checkCool(file = coolFile), ".cool")
    
    # Resolutions
    expect_equal(readCoolBpResolutions(file = coolFile), 1000000)
    expect_error(readCoolNorms(file = coolFile, resolution = 5000), "resolution")
    expect_error(readCoolChroms(file = coolFile, resolution = 5000), "resolution")
    
    # Norms
    expect_equal(readCoolNorms(file = coolFile), c("BALANCE", "NONE"))
    
    # Chroms
    cool_chroms <- readCoolChroms(file = coolFile) |>
        pull(name)
    expect_setequal(cool_chroms, paste0("chr", c(seq(1:22), "X", "Y", "M")))
    expect_error(.checkCoolErrors(file = coolFile, chrom = "2",
                                  chromstart = NULL, chromend = NULL,
                                  zrange = NULL, altchrom = NULL,
                                  altchromstart = NULL, altchromend = NULL,
                                  norm = "NONE", resolution = 1000000), 
                 "chrom")
    
    ## .mcool file
    mcoolFile <- file.path(tempdir(), "LEUK_HEK_PJA27_inter_30.mcool")
    download.file(url = "https://zenodo.org/records/10906240/files/LEUK_HEK_PJA27_inter_30.mcool?download=1",
                  destfile = mcoolFile, mode = "wb")
    on.exit(unlink(mcoolFile))
    
    # File type
    expect_equal(.checkCool(file = mcoolFile), ".mcool")
    
    # Resolutions
    expect_setequal(readCoolBpResolutions(file = mcoolFile),
                    c(5000, 10000, 25000, 50000, 100000, 250000, 500000,
                      1000000, 2500000))
    expect_error(readCoolNorms(file = mcoolFile, resolution = 2500), 
                 "resolution")
    
    # Norms
    # Vector of norms for one resolution
    expect_setequal(readCoolNorms(file = mcoolFile, resolution = 10000),
                    c("KR", "VC", "VC_SQRT", "NONE"))
    # List of vectors for all resolutions
    expect_type(readCoolNorms(file = mcoolFile), "list")
    
    # Chroms
    expect_type(readCoolChroms(file = mcoolFile), "list")
    expect_error(.checkCoolErrors(file = mcoolFile, chrom = "chr2",
                                  chromstart = NULL, chromend = NULL,
                                  zrange = NULL, altchrom = NULL,
                                  altchromstart = NULL, altchromend = NULL,
                                  norm = "NONE", resolution = 1000000), 
                 "chrom")
    
})
