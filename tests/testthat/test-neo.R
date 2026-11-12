test_that("readNeoAssembly parses fragments, offsets, and flips", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )

    neo <- readNeoAssembly(assembly, "C1", extend = 0, quiet = TRUE)

    expect_s3_class(neo, "neoAssembly")
    expect_equal(neo$total_length, 800)
    expect_equal(nrow(neo$fragments), 2)
    expect_equal(neo$fragments$chrom, c("chr1", "chr2"))
    expect_equal(neo$fragments$flipped, c(TRUE, FALSE))
    expect_equal(neo$fragments$offset, c(0, 400))
    expect_equal(neo$fragments$length, c(400, 400))
})

test_that("readNeoAssembly extends only the first and last fragments", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )

    neo <- readNeoAssembly(assembly, "C1", extend = 100, quiet = TRUE)

    expect_equal(neo$fragments$length, c(500, 500))
    expect_equal(neo$fragments$offset, c(0, 500))
    expect_equal(neo$total_length, 1000)
})

test_that("mapNeoFeatures flips coordinates and strand", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    features <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(0, 0),
        end = c(100, 100),
        name = c("f1", "f2"),
        strand = c("+", "+")
    )
    mapped <- mapNeoFeatures(features, neo)

    expect_equal(mapped$chrom, c("neo", "neo"))
    expect_equal(mapped$fragment, c(1L, 2L))
    expect_equal(mapped$start, c(300, 400))
    expect_equal(mapped$end, c(400, 500))
    expect_equal(mapped$strand, c("-", "+"))
})

test_that("readNeoHic reconstructs a mcool matrix in neo coordinates", {
    mcool <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(mcool)

    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    neoHic <- readNeoHic(
        mcool, neo, resolution = 100, norm = "NONE", quiet = TRUE
    )

    expect_s3_class(neoHic, "neoHic")
    expect_equal(neoHic$resolution, 100)
    expect_equal(neoHic$neo$total_length, 800)
    expect_equal(nrow(neoHic$pair_counts), 3)
    expect_true(all(c("fragment1", "fragment2", "pixels", "counts") %in%
                    colnames(neoHic$pair_counts)))
    expect_true(all(neoHic$matrix$x <= neoHic$matrix$y))
    expect_false(any(duplicated(neoHic$matrix[, c("x", "y")])))

    expected <- fx$pixels[
        (fx$pixels$bin1_id %in% 0:3 & fx$pixels$bin2_id %in% 0:3) |
        (fx$pixels$bin1_id %in% 0:3 & fx$pixels$bin2_id %in% 10:13) |
        (fx$pixels$bin1_id %in% 10:13 & fx$pixels$bin2_id %in% 10:13),
    ]
    expect_equal(sum(neoHic$matrix$counts), sum(expected$count))

    getCount <- function(x, y) {
        neoHic$matrix$counts[neoHic$matrix$x == x & neoHic$matrix$y == y]
    }
    ## chr1 bin 0 (flipped to neo 300) vs chr2 bin 10 (neo 400)
    expect_equal(
        getCount(300, 400),
        fx$pixels$count[fx$pixels$bin1_id == 0 & fx$pixels$bin2_id == 10]
    )
    ## chr1 bin 1 (flipped to neo 200) vs chr2 bin 11 (neo 500)
    expect_equal(
        getCount(200, 500),
        fx$pixels$count[fx$pixels$bin1_id == 1 & fx$pixels$bin2_id == 11]
    )
})

test_that("plotNeoHic creates a plot object", {
    mcool <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(mcool)

    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    neoHic <- readNeoHic(
        mcool, neo, resolution = 100, norm = "NONE", quiet = TRUE
    )

    pageCreate(width = 4, height = 3, default.units = "inches",
               showGuides = FALSE)
    plot <- plotNeoHic(
        neoHic, type = "triangle",
        x = 0.5, y = 0.5, width = 2, height = 1.5,
        just = c("left", "top"), draw = FALSE
    )
    expect_s3_class(plot, "hicTriangle")
})

test_that("readNeoHic accepts an assembly path and ID", {
    mcool <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(mcool)

    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )

    neoHic <- readNeoHic(
        mcool, assembly, id = "C1",
        resolution = 100, norm = "NONE", quiet = TRUE
    )
    expect_s3_class(neoHic, "neoHic")
    expect_equal(neoHic$neo$id, "C1")
    expect_true(nrow(neoHic$matrix) > 0)
})

test_that("readNeoHic accepts precomputed genomic pair data", {
    mcool <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(mcool)

    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    binChrom <- ifelse(fx$pixels$bin1_id < fx$nBins1, "chr1", "chr2")
    binChrom2 <- ifelse(fx$pixels$bin2_id < fx$nBins1, "chr1", "chr2")
    pairData <- data.frame(
        chrom1 = binChrom,
        start1 = fx$starts[fx$pixels$bin1_id + 1],
        chrom2 = binChrom2,
        start2 = fx$starts[fx$pixels$bin2_id + 1],
        counts = fx$pixels$count
    )
    pairData <- pairData[
        (pairData$chrom1 == "chr1" & pairData$start1 < 400 &
         pairData$chrom2 == "chr1" & pairData$start2 < 400) |
        (pairData$chrom1 == "chr1" & pairData$start1 < 400 &
         pairData$chrom2 == "chr2" & pairData$start2 < 400) |
        (pairData$chrom1 == "chr2" & pairData$start1 < 400 &
         pairData$chrom2 == "chr2" & pairData$start2 < 400),
    ]

    fromFile <- readNeoHic(
        mcool, neo, resolution = 100, norm = "NONE", quiet = TRUE
    )
    fromDf <- readNeoHic(
        pairData, neo, resolution = 100, quiet = TRUE
    )

    expect_equal(fromDf$matrix, fromFile$matrix)
})

test_that("annoNeoFragmentLabels supports coordinate and inversion labels", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    labels <- annoNeoFragmentLabels(
        neo,
        x = 0.5, y = 0.5, width = 3,
        label.style = "coordinate",
        label.position = "below",
        label.offset = 0.05,
        label.fontface = "italic",
        show.inversion = TRUE
    )
    expect_type(labels, "list")
    expect_length(labels, 2)
})

test_that("plotNeoHic can flip the x-axis horizontally", {
    mcool <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(mcool)
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    neoHic <- readNeoHic(mcool, neo, resolution = 100,
                          norm = "NONE", quiet = TRUE)

    pageCreate(width = 3, height = 3, default.units = "inches",
               showGuides = FALSE)
    hic <- plotNeoHic(
        neoHic, type = "square", flip = "horizontal",
        x = 0.5, y = 0.5, width = 2, height = 2,
        draw = FALSE, zrange = c(0, 2)
    )
    expect_equal(hic$xaxisRange, c(neo$total_length, 0))
    expect_equal(hic$xaxisDir, "decrease")

    tri <- plotNeoHic(
        neoHic, type = "triangle", flip = "horizontal",
        x = 0.5, y = 0.5, width = 2, height = 1,
        draw = FALSE, zrange = c(0, 2)
    )
    expect_equal(tri$xaxisRange, c(neo$total_length, 0))
    expect_equal(tri$xaxisDir, "decrease")
})

test_that("plotNeoGenes draws IGV-style gene bodies and exon blocks", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    genes <- data.frame(
        chrom = c("chr1", "chr1", "chr2", "chr2"),
        start = c(100, 150, 100, 300),
        end = c(120, 200, 180, 350),
        strand = c("+", "+", "-", "-"),
        gene = c("GENE1", "GENE1", "GENE2", "GENE3"),
        feature = "exon",
        gene_type = c("protein_coding", "protein_coding",
                      "protein_coding", "lncRNA")
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 3, height = 0.5,
        labels = "GENE1", geneHighlights = "GENE1",
        gene.type = "protein_coding", show.strand = TRUE, max.rows = 2
    )

    ## lncRNA GENE3 is filtered out by gene.type
    expect_setequal(unique(res$genes$gene), c("GENE1", "GENE2"))
    ## Each gene is one row with a list of reduced exon intervals
    expect_equal(nrow(res$genes), 2)
    expect_equal(length(res$genes$exonStart[[1]]), 2)
    ## IGV style: body segments plus exon rectangles
    classes <- vapply(res$grobs, function(g) class(g)[1], character(1))
    expect_true("segments" %in% classes)
    expect_true("rect" %in% classes)
    ## Highlighted gene is red, non-highlighted gene stays black
    expect_equal(res$grobs[[1]]$gp$col, "red")
    grobCols <- vapply(res$grobs, function(g) {
        if (is.null(g$gp$col)) NA_character_ else as.character(g$gp$col)[1]
    }, character(1))
    expect_true("red" %in% grobCols)
    expect_true("black" %in% grobCols)

    ## Exon blocks are centred on the gene body, not offset above it
    body <- res$grobs[[1]]
    exon <- res$grobs[[which(classes == "rect")[1]]]
    expect_equal(as.numeric(exon$y),
                 as.numeric(body$y0) + as.numeric(exon$height) / 2)
    expect_equal(exon$just, c(0, 1))
})

test_that("plotNeoGenes can hide strand labels", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = 100, end = 200, strand = "+",
        gene = "GENE1", feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    noStrand <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 3, height = 0.3,
        labels = FALSE, show.strand = FALSE
    )
    strand <- plotNeoGenes(
        genes, neo, x = 0.5, y = 1, width = 3, height = 0.3,
        labels = FALSE, show.strand = TRUE
    )
    expect_false(any(vapply(noStrand$grobs,
                            function(g) identical(g$label, "+"), logical(1))))
    expect_true(any(vapply(strand$grobs,
                           function(g) identical(g$label, "+"), logical(1))))
})

test_that("plotNeoAssembly and annoNeoDomains expose seam styling", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    pageCreate(width = 4, height = 3, default.units = "inches",
               showGuides = FALSE)
    asm <- plotNeoAssembly(
        neo, x = 0.5, y = 0.5, width = 3, height = 0.2,
        labels = FALSE, seam.color = "purple", seam.lwd = 2,
        seam.lty = 3
    )
    expect_equal(asm$seams[[1]]$gp$col, "purple")
    expect_equal(asm$seams[[1]]$gp$lwd, 2)
    expect_equal(asm$seams[[1]]$gp$lty, 3)

    domainPlot <- plotNeoHic(
        readNeoHic(
            {
                mcool <- withr::local_tempfile(fileext = ".mcool")
                make_test_mcool(mcool)
                mcool
            },
            neo, resolution = 100, norm = "NONE", quiet = TRUE
        ),
        type = "triangle",
        x = 0.5, y = 1, width = 3, height = 1,
        draw = FALSE, zrange = c(0, 2)
    )
    domains <- annoNeoDomains(
        domainPlot, neo, seam.color = "green", seam.lwd = 1.5,
        seam.lty = 2
    )
    expect_equal(length(domains$seams), 1)
    seamKids <- domains$seams[[1]]$grobs$children
    expect_true(length(seamKids) > 0)
    expect_equal(seamKids[[1]]$gp$col, "green")
    expect_equal(seamKids[[1]]$gp$lwd, 1.5)
    expect_equal(seamKids[[1]]$gp$lty, 2)
})

test_that("annoNeoFragmentLabels accepts explicit per-fragment colors", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    labels <- annoNeoFragmentLabels(
        neo, x = 0.5, y = 0.5, width = 3,
        labels = c("frag1", "frag2"),
        label.colors = c("red", "green")
    )
    expect_equal(labels[[1]]$gp$col, "red")
    expect_equal(labels[[2]]$gp$col, "green")
})

test_that("plotNeoGenes reads arbitrary GTF files", {
    skip_if_not_installed("rtracklayer")
    gtf <- withr::local_tempfile(fileext = ".gtf")
    writeLines(c(
        'chr1\tTEST\tgene\t100\t200\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr1\tTEST\texon\t100\t120\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr1\tTEST\texon\t180\t200\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr1\tTEST\tgene\t300\t400\t.\t+\t.\tgene_id "G2"; gene_name "GENE2"; gene_type "lncRNA";',
        'chr1\tTEST\texon\t300\t400\t.\t+\t.\tgene_id "G2"; gene_name "GENE2"; gene_type "lncRNA";'
    ), gtf)
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        gtf, neo, x = 0.5, y = 0.5, width = 3, height = 0.3,
        labels = "GENE1", gene.type = "protein_coding"
    )
    expect_equal(unique(res$genes$gene), "GENE1")
    expect_equal(length(res$genes$exonStart[[1]]), 2)
})

test_that("collapsed display mode draws solid gene blocks without exon ticks", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 300), end = c(200, 400),
        strand = "+", gene = c("GENE1", "GENE2"), feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.4,
        displayMode = "collapsed", labels = FALSE,
        show.strand = "none", max.rows = 2
    )

    expect_equal(res$displayMode, "collapsed")
    classes <- vapply(res$grobs, function(g) class(g)[1], character(1))
    expect_true(all(classes == "rect"))
    minPx <- convertWidth(unit(3, "points"), unitTo = "inches",
                          valueOnly = TRUE)
    widths <- vapply(res$grobs, function(g) as.numeric(g$width), numeric(1))
    expect_true(all(widths >= minPx - 1e-8))
})

test_that("expanded display mode draws body lines and exon blocks", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 150), end = c(120, 220),
        strand = "+", gene = "GENE1", feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.4,
        displayMode = "expanded", labels = FALSE,
        show.strand = "none", max.rows = 1
    )

    expect_equal(res$displayMode, "expanded")
    classes <- vapply(res$grobs, function(g) class(g)[1], character(1))
    expect_true("segments" %in% classes)
    expect_true("rect" %in% classes)
    exon <- res$grobs[[which(classes == "rect")[1]]]
    expect_equal(exon$just, c(0, 1))
    minPx <- convertWidth(unit(1, "points"), unitTo = "inches",
                          valueOnly = TRUE)
    expect_gte(as.numeric(exon$width), minPx - 1e-8)
})

test_that("auto display mode switches on basepairs per inch", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = 100, end = 200, strand = "+",
        gene = "GENE1", feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    collapsed <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "auto", auto.collapse.bp.per.inch = 1,
        labels = FALSE, show.strand = "none"
    )
    expanded <- plotNeoGenes(
        genes, neo, x = 0.5, y = 1, width = 2, height = 0.3,
        displayMode = "auto", auto.collapse.bp.per.inch = 1e9,
        labels = FALSE, show.strand = "none"
    )
    expect_equal(collapsed$displayMode, "collapsed")
    expect_equal(expanded$displayMode, "expanded")
})

test_that("chromstart and chromend clip and rescale the neo window", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 300), end = c(200, 400),
        strand = "+", gene = c("GENE1", "GENE2"), feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "collapsed", chromstart = 620, chromend = 800,
        labels = FALSE, show.strand = "none"
    )
    expect_equal(res$genes$gene, "GENE2")
    expect_equal(res$chromstart, 620)
    expect_equal(res$chromend, 800)
    ## chr2 fragment starts at neo 400, so 300 bp is neo 700. The window
    ## starts at 620 and spans 180 bp.
    expect_equal(as.numeric(res$grobs[[1]]$x), 0.5 + 80 / 180 * 2)
})

test_that("plotNeoGenes warns and drops genes when rows overflow", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 100, 100), end = c(300, 300, 300),
        strand = "+", gene = c("G1", "G2", "G3"), feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    expect_warning(
        res <- plotNeoGenes(
            genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
            displayMode = "collapsed", labels = FALSE,
            show.strand = "none", max.rows = 1
        ),
        "Not enough rows"
    )
    expect_equal(nrow(res$genes), 1)
})

test_that("plotNeoGenes suppresses overlapping labels", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 150), end = c(140, 190),
        strand = "+", gene = c("GENE1", "GENE2"), feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "collapsed", labels = c("GENE1", "GENE2"),
        show.strand = "none", label.overlap = TRUE, max.rows = 1
    )
    texts <- Filter(function(g) class(g)[1] == "text", res$grobs)
    labels <- vapply(texts, function(g) as.character(g$label), character(1))
    ys <- vapply(texts, function(g) as.numeric(g$y), numeric(1))
    expect_equal(sort(labels), c("GENE1", "GENE2"))
    expect_false(isTRUE(all.equal(ys[1], ys[2])))
})

test_that("single-exon genes are solid blocks in expanded mode", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = 100, end = 200, strand = "+",
        gene = "GENE1", feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", labels = FALSE, show.strand = "none"
    )
    classes <- vapply(res$grobs, function(g) class(g)[1], character(1))
    expect_false("segments" %in% classes)
    expect_true("rect" %in% classes)
})

test_that("UTR and CDS blocks use half and full heights", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2",
        start = c(100, 160, 120, 160),
        end = c(140, 220, 140, 200),
        strand = "+",
        gene = "GENE1",
        feature = c("exon", "exon", "CDS", "CDS")
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", labels = FALSE, show.strand = "none",
        max.rows = 1
    )
    classes <- vapply(res$grobs, function(g) class(g)[1], character(1))
    heights <- vapply(res$grobs[classes == "rect"],
                      function(g) as.numeric(g$height), numeric(1))
    expect_equal(length(heights), 4)
    expect_equal(heights[1] / heights[3], 0.5)
    expect_equal(heights[2] / heights[4], 0.5)
})

test_that("strand arrows are emitted in collapsed mode", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = 100, end = 300, strand = "+",
        gene = "GENE1", feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "collapsed", labels = FALSE,
        show.strand = "arrows"
    )
    classes <- vapply(res$grobs, function(g) class(g)[1], character(1))
    expect_true("polygon" %in% classes)
})

test_that("plotNeoGenes reads CDS and UTR features from a GTF", {
    skip_if_not_installed("rtracklayer")
    gtf <- withr::local_tempfile(fileext = ".gtf")
    writeLines(c(
        'chr2\tTEST\tgene\t100\t300\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\texon\t100\t140\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\texon\t200\t300\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\tCDS\t120\t140\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\tCDS\t200\t280\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\tUTR\t100\t120\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\tUTR\t280\t300\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";'
    ), gtf)
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        gtf, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", labels = FALSE, show.strand = "none"
    )
    expect_equal(nrow(res$genes), 1)
    expect_equal(length(res$genes$cdsStart[[1]]), 2)
    expect_equal(length(res$genes$exonStart[[1]]), 2)
})

test_that("plotNeoGenes accepts GRanges input", {
    skip_if_not_installed("GenomicRanges")
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    gr <- GenomicRanges::GRanges(
        "chr2",
        IRanges::IRanges(100, 200),
        strand = "+",
        gene = "GENE1",
        feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        gr, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", labels = FALSE, show.strand = "none"
    )
    expect_equal(nrow(res$genes), 1)
    expect_equal(res$genes$gene, "GENE1")
})

test_that("squished display mode uses shorter rows than expanded", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = 100, end = 200, strand = "+",
        gene = "GENE1", feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    expanded <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.4,
        displayMode = "expanded", labels = FALSE, show.strand = "none"
    )
    squished <- plotNeoGenes(
        genes, neo, x = 0.5, y = 1, width = 2, height = 0.4,
        displayMode = "squished", labels = FALSE, show.strand = "none"
    )
    expect_equal(squished$displayMode, "squished")
    expandedRect <- expanded$grobs[[which(vapply(
        expanded$grobs, function(g) class(g)[1] == "rect", logical(1)
    ))[1]]]
    squishedRect <- squished$grobs[[which(vapply(
        squished$grobs, function(g) class(g)[1] == "rect", logical(1)
    ))[1]]]
    expect_lt(as.numeric(squishedRect$height), as.numeric(expandedRect$height))
})

test_that("geneBackground controls non-highlighted genes", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines(
        "C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400",
        assembly
    )
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 300), end = c(200, 400),
        strand = "+", gene = c("GENE1", "GENE2"), feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "collapsed", labels = FALSE, show.strand = "none",
        geneHighlights = "GENE1", geneBackground = "grey"
    )
    cols <- vapply(res$grobs, function(g) {
        if (is.null(g$gp$col)) NA_character_ else as.character(g$gp$col)[1]
    }, character(1))
    expect_true("red" %in% cols)
    expect_true("grey" %in% cols)
})

test_that("geneOrder prioritizes labels", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 150), end = c(140, 190),
        strand = "+", gene = c("LONGNAME1", "LONGNAME2"), feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "collapsed", labels = c("LONGNAME1", "LONGNAME2"),
        geneOrder = c("LONGNAME2", "LONGNAME1"), show.strand = "none",
        label.overlap = TRUE, max.rows = 1
    )
    texts <- Filter(function(g) class(g)[1] == "text", res$grobs)
    labels <- vapply(texts, function(g) as.character(g$label), character(1))
    ys <- vapply(texts, function(g) as.numeric(g$y), numeric(1))
    ## Priority label stays on the lower level; the other is stacked above.
    ## plotText stores bottom-based y, so the priority label sits lower.
    expect_lt(ys[match("LONGNAME2", labels)],
              ys[match("LONGNAME1", labels)])
})

test_that("labelClusters draws rotated labels without stacking", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 150), end = c(140, 190),
        strand = "+", gene = c("LONGNAME1", "LONGNAME2"), feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "collapsed", labels = c("LONGNAME1", "LONGNAME2"),
        labelClusters = TRUE, labelClusterRot = 45,
        show.strand = "none", label.overlap = TRUE, max.rows = 1
    )
    texts <- Filter(function(g) class(g)[1] == "text", res$grobs)
    expect_equal(length(texts), 2)
    ys <- vapply(texts, function(g) as.numeric(g$y), numeric(1))
    rots <- vapply(texts, function(g) as.numeric(g$grobs$rot), numeric(1))
    expect_equal(length(unique(ys)), 1)
    expect_true(all(rots == 45))
})

test_that("plotNeoGenes supports TxDb and OrgDb annotation", {
    skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    skip_if_not_installed("org.Hs.eg.db")
    suppressPackageStartupMessages({
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        library(org.Hs.eg.db)
    })
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\t12,65820000\t12,65830000", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        data = TxDb.Hsapiens.UCSC.hg38.knownGene,
        orgdb = org.Hs.eg.db,
        neo = neo, x = 0.5, y = 0.5, width = 3, height = 0.5,
        displayMode = "expanded", labels = "HMGA2", show.strand = "none"
    )
    expect_true(any(res$genes$gene == "HMGA2"))

    res2 <- plotNeoGenes(
        txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
        orgdb = org.Hs.eg.db,
        neo = neo, x = 0.5, y = 1.1, width = 3, height = 0.5,
        displayMode = "collapsed", labels = FALSE, show.strand = "none"
    )
    expect_true(any(res2$genes$gene == "HMGA2"))
})

test_that("GTF readers agree and cache the parsed annotation", {
    skip_if_not_installed("rtracklayer")
    gtf <- withr::local_tempfile(fileext = ".gtf")
    writeLines(c(
        'chr2\tTEST\tgene\t100\t300\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\texon\t100\t140\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\texon\t200\t300\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";'
    ), gtf)
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    fast <- plotNeoGenes(
        gtf, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", labels = FALSE, show.strand = "none",
        reader = "readGFF"
    )
    fallback <- plotNeoGenes(
        gtf, neo, x = 0.5, y = 1, width = 2, height = 0.3,
        displayMode = "expanded", labels = FALSE, show.strand = "none",
        reader = "import"
    )
    expect_equal(fast$genes$gene, fallback$genes$gene)
    expect_equal(fast$genes$start, fallback$genes$start)
    expect_equal(fast$genes$end, fallback$genes$end)

    cache <- withr::local_tempdir()
    plotNeoGenes(
        gtf, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "collapsed", labels = FALSE, show.strand = "none",
        reader = "readGFF", cache = cache
    )
    expect_equal(length(list.files(cache)), 1)
})

test_that("tabix reader returns the same annotation", {
    skip_if_not_installed("rtracklayer")
    skip_if_not_installed("Rsamtools")
    gtf <- withr::local_tempfile(fileext = ".gtf")
    writeLines(c(
        'chr2\tTEST\tgene\t100\t300\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\texon\t100\t140\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";',
        'chr2\tTEST\texon\t200\t300\t.\t+\t.\tgene_id "G1"; gene_name "GENE1"; gene_type "protein_coding";'
    ), gtf)
    ## Tabix requires position-sorted input.
    lines <- readLines(gtf)
    fields <- strsplit(lines, "\t", fixed = TRUE)
    ord <- order(
        vapply(fields, `[`, character(1), 1),
        as.numeric(vapply(fields, `[`, character(1), 4))
    )
    sorted <- withr::local_tempfile(fileext = ".gtf")
    writeLines(lines[ord], sorted)
    gz <- paste0(sorted, ".gz")
    withr::defer(unlink(c(gz, paste0(gz, ".tbi"))))
    Rsamtools::bgzip(sorted, dest = gz, overwrite = TRUE)
    Rsamtools::indexTabix(gz, format = "gff")

    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        gz, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", labels = FALSE, show.strand = "none",
        reader = "tabix"
    )
    expect_equal(res$genes$gene, "GENE1")
    expect_equal(length(res$genes$exonStart[[1]]), 2)
})

test_that("focus mode fades non-priority genes but keeps their exons", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2",
        start = c(100, 220, 320, 370),
        end = c(140, 300, 360, 390),
        strand = "+",
        gene = c("GENE1", "GENE1", "GENE2", "GENE2"),
        feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", geneMode = "focus", focus.alpha = 0.5,
        geneHighlights = "GENE1", geneBackground = "grey",
        labels = "GENE1", show.strand = "arrows", max.rows = 1
    )
    expect_equal(res$geneMode, "focus")
    rects <- Filter(function(g) class(g)[1] == "rect", res$grobs)
    cols <- vapply(rects, function(g) as.character(g$gp$col)[1], character(1))
    alphas <- vapply(rects, function(g) as.numeric(g$gp$alpha)[1], numeric(1))
    expect_true(any(cols == "red" & alphas == 1))
    expect_true(any(cols == "grey" & alphas == 0.5))
    segments <- Filter(function(g) class(g)[1] == "segments", res$grobs)
    segCols <- vapply(segments, function(g) as.character(g$gp$col)[1],
                      character(1))
    segAlphas <- vapply(segments, function(g) as.numeric(g$gp$alpha)[1],
                        numeric(1))
    expect_true(any(segCols == "grey" & segAlphas == 0.5))
    texts <- Filter(function(g) class(g)[1] == "text", res$grobs)
    expect_equal(length(texts), 1)
    expect_equal(texts[[1]]$label, "GENE1")
})

test_that("compact mode draws all genes in a single row", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2",
        start = c(100, 220, 320, 370),
        end = c(140, 300, 360, 390),
        strand = "+",
        gene = c("GENE1", "GENE1", "GENE2", "GENE2"),
        feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", geneMode = "compact",
        geneHighlights = "GENE1", show.strand = "none", max.rows = 4
    )
    expect_equal(res$geneMode, "compact")
    expect_true(all(res$genes$row == 1))
    rects <- Filter(function(g) class(g)[1] == "rect", res$grobs)
    cols <- vapply(rects, function(g) as.character(g$gp$col)[1], character(1))
    expect_true("red" %in% cols)
    expect_true("grey" %in% cols)
})

test_that("focus mode keeps overflowing background genes on the last row", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,+,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    genes <- data.frame(
        chrom = "chr2", start = c(100, 100, 100), end = c(300, 300, 300),
        strand = "+", gene = c("GENE1", "GENE2", "GENE3"), feature = "exon"
    )

    pageCreate(width = 4, height = 2, default.units = "inches",
               showGuides = FALSE)
    res <- plotNeoGenes(
        genes, neo, x = 0.5, y = 0.5, width = 2, height = 0.3,
        displayMode = "expanded", geneMode = "focus",
        geneHighlights = "GENE1", labels = FALSE,
        show.strand = "none", max.rows = 1
    )
    expect_equal(nrow(res$genes), 3)
    expect_true(all(res$genes$row == 1))
})

test_that("getNeoSignal maps and bins bigWig-like signal", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    signal <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100, 100), end = c(200, 200),
        score = c(2, 4)
    )
    res <- getNeoSignal(signal, neo, binSize = 100)
    expect_equal(res$chrom, rep("neo", 2))
    expect_equal(res$start, c(200, 500))
    expect_equal(res$end, c(299, 599))
    expect_equal(res$score, c(2, 4))

    gr <- GenomicRanges::GRanges(
        c("chr1", "chr2"),
        IRanges::IRanges(c(100, 100), c(200, 200)),
        score = c(2, 4)
    )
    grRes <- getNeoSignal(gr, neo, binSize = 100)
    expect_equal(grRes$start, c(200, 500))

    window <- getNeoSignal(signal, neo, binSize = 100,
                           chromstart = 450, chromend = 650)
    expect_equal(window$start, 500)
})

test_that("getNeoV4C maps virtual-4C tables", {
    assembly <- withr::local_tempfile(fileext = ".assemblies.txt")
    writeLines("C1\ttranslocation,1,0,-,2,0,+\t1,400\t2,400", assembly)
    neo <- readNeoAssembly(assembly, "C1", quiet = TRUE)
    tsv <- withr::local_tempfile(fileext = ".tsv")
    writeLines(c(
        "chrom\tstart\tend\tbalanced",
        "chr2\t100\t200\t1.5",
        "chr2\t200\t300\t2.5"
    ), tsv)
    res <- getNeoV4C(tsv, neo)
    expect_equal(res$chrom, rep("neo", 2))
    expect_equal(res$start, c(500, 600))
    expect_equal(res$score, c(1.5, 2.5))
    expect_true(all(res$end < c(600, 700)))
})
