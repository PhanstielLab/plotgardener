## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    fig.align = "center",
    fig.width = 3,
    fig.height = 3,
    collapse = TRUE,
    comment = "#>",
    fig.retina = 1,
    warning = FALSE,
    message = FALSE
)
library(grid)
library(BentoBox)

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  wholeFile <- bb_readBigwig("/path/to/bigWig")
#  
#  region <- bb_readBigwig("/path/to/bigWig",
#      chrom = "chr1",
#      chromstart = 1000000, chromend = 2000000
#  )
#  
#  regionPlus <- bb_readBigwig("/path/to/bigWig",
#      chrom = "chr1",
#      chromstart = 1000000, chromend = 2000000,
#      strand = "+"
#  )

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  chrom <- bb_readHic("/path/to/hic",
#      chrom = "chr1",
#      resolution = 250000, res_scale = "BP", norm = "NONE"
#  )
#  
#  chromRegion <- bb_readHic("/path/to/hic",
#      chrom = "chr1",
#      chromstart = 1000000, chromend = 2000000,
#      resolution = 10000, res_scale = "BP", norm = "KR"
#  )
#  
#  twoChroms <- bb_readHic("/path/to/hic",
#      chrom = "chr1", altchrom = "chr2",
#      resolution = 250000, res_scale = "BP"
#  )

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  library(data.table)
#  data <- data.table::fread("/path/to/file")
#  
#  library(rtracklayer)
#  data <- rtracklayer::import(con = "/path/to/file", format = "fileFormat")

## ----hic_quickplot, eval=TRUE, echo=TRUE, message=FALSE-----------------------
## Load BentoBox
library(BentoBox)

## Load example Hi-C data
data("bb_imrHicData")

## Quick plot Hi-C data
bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000
)

## ----signal_quickplot, eval=TRUE, echo=TRUE, message=FALSE--------------------
## Load BentoBox
library(BentoBox)

## Load example signal data
data("bb_imrH3K27acData")

## Quick plot signal data
bb_plotSignal(
    data = bb_imrH3K27acData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000
)

## ----gene_quickplot, eval=TRUE, echo=TRUE, message=FALSE----------------------
## Load BentoBox
library(BentoBox)

## Load hg19 genomic annotation packages
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

## Quick plot genes
bb_plotGenes(
    assembly = "hg19",
    chrom = "chr21", chromstart = 28000000, chromend = 30300000
)

## ----gwas_quickplot, eval=TRUE, echo=TRUE, message=FALSE----------------------
## Load BentoBox
library(BentoBox)

## Load hg19 genomic annotation packages
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## Load example GWAS data
data("bb_gwasData")

## Quick plot GWAS data
bb_plotManhattan(
    data = bb_gwasData, fill = c("steel blue", "grey"),
    ymax = 1.1, cex = 0.20
)

## ----quickpage, echo=TRUE, fig.height=4, fig.width=4, message=FALSE-----------
bb_pageCreate(width = 3.25, height = 3.25, default.units = "inches")

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  data("bb_imrHicData")
#  hicPlot <- bb_plotHicSquare(
#      data = bb_imrHicData,
#      chrom = "chr21", chromstart = 28000000, chromend = 30300000,
#      x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
#  )

## ----quickpageHic, echo=FALSE, fig.height=4, fig.width=4, message=FALSE-------
bb_pageCreate(width = 3.25, height = 3.25, default.units = "inches")
data("bb_imrHicData")
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
)

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  bb_annoHeatmapLegend(
#      plot = hicPlot,
#      x = 2.85, y = 0.25, width = 0.1, height = 1.25, default.units = "inches"
#  )
#  
#  bb_annoGenomeLabel(
#      plot = hicPlot,
#      x = 0.25, y = 2.75, width = 2.5, height = 0.25, default.units = "inches"
#  )

## ----quickpageAnno, echo=FALSE, fig.height=4, fig.width=4, message=FALSE------
bb_pageCreate(width = 3.25, height = 3.25, default.units = "inches")
data("bb_imrHicData")
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
)
bb_annoHeatmapLegend(
    plot = hicPlot,
    x = 2.85, y = 0.25, width = 0.1, height = 1.25, default.units = "inches"
)

bb_annoGenomeLabel(
    plot = hicPlot,
    x = 0.25, y = 2.75, width = 2.5, height = 0.25, default.units = "inches"
)

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  bb_pageGuideHide()

## ----quickpageHide, echo=FALSE, fig.height=4, fig.width=4, message=FALSE------
bb_pageCreate(
    width = 3.25, height = 3.25, default.units = "inches",
    xgrid = 0, ygrid = 0, showGuides = FALSE
)
data("bb_imrHicData")
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
)
bb_annoHeatmapLegend(
    plot = hicPlot,
    x = 2.85, y = 0.25, width = 0.1, height = 1.25, default.units = "inches"
)

bb_annoGenomeLabel(
    plot = hicPlot,
    x = 0.25, y = 2.75, width = 2.5, height = 0.25, default.units = "inches"
)

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  pdf(width = 3.25, height = 3.25)
#  
#  bb_pageCreate(width = 3.25, height = 3.25, default.units = "inches")
#  data("bb_imrHicData")
#  hicPlot <- bb_plotHicSquare(
#      data = bb_imrHicData,
#      chrom = "chr21", chromstart = 28000000, chromend = 30300000,
#      x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
#  )
#  bb_annoHeatmapLegend(
#      plot = hicPlot,
#      x = 2.85, y = 0.25, width = 0.1, height = 1.25, default.units = "inches"
#  )
#  
#  bb_annoGenomeLabel(
#      plot = hicPlot,
#      x = 0.25, y = 2.75, width = 2.5, height = 0.25, default.units = "inches"
#  )
#  bb_pageGuideHide()
#  
#  dev.off()

## ----page_demo01, fig.width=11, fig.height=12---------------------------------
bb_pageCreate(width = 8.5, height = 11, default.units = "inches")

## ----page_demo02, fig.width=5, fig.height=4-----------------------------------
bb_pageCreate(width = 8, height = 8, xgrid = 1, ygrid = 1, default.units = "cm")

## ----page_demo03, fig.width=5, fig.height=4-----------------------------------
bb_pageCreate(
    width = 3, height = 3, xgrid = 0, ygrid = 0,
    default.units = "inches"
)

## ----page_demo04, fig.width=5, fig.height=4-----------------------------------
bb_pageCreate(width = 3, height = 3, default.units = "inches")
## Add a horizontal guide at y = 2.25 inches
bb_pageGuideHorizontal(y = 2.25, default.units = "inches")
## Add a vertical guide at x = 0.75 inches
bb_pageGuideVertical(x = 0.75, default.units = "inches")

## ----page_demo05, fig.width=5, fig.height=4-----------------------------------
## Create page
bb_pageCreate(width = 3, height = 3, default.units = "inches")
## Remove guides
bb_pageGuideHide()

## ----coord_01, eval=TRUE, message=FALSE,echo=TRUE, fig.width=5, fig.height=4----
bb_pageCreate(width = 3, height = 3, default.units = "inches")

## ----eval=FALSE, message=FALSE,echo=TRUE--------------------------------------
#  bb_plotRect(
#      x = unit(0.5, "npc"), y = unit(0.5, "npc"), width = 1, height = 1,
#      default.units = "inches"
#  )

## ----coord_02, eval=TRUE, message=FALSE,echo=FALSE, fig.width=5, fig.height=4----
bb_pageCreate(width = 3, height = 3, default.units = "inches")
bb_plotRect(x = unit(0.5, "npc"), y = unit(0.5, "npc"), width = 1, height = 1)

## ----coord_03, echo=TRUE, fig.height=2.5, fig.width=7, message=FALSE----------
bb_pageCreate(
    width = 5, height = 1.5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

data("bb_imrH3K27acData")
signalPlot <- bb_plotSignal(
    data = bb_imrH3K27acData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 0.25, width = 4, height = 0.75, default.units = "inches"
)
bb_annoGenomeLabel(plot = signalPlot, x = 0.5, y = 1.01)

## Annotate text at average x-coordinate of data peak
peakScore <- bb_imrH3K27acData[which(
    bb_imrH3K27acData$score == max(bb_imrH3K27acData$score)
), ]
peakPos <- round((min(peakScore$start) + max(peakScore$end)) * 0.5)
bb_annoText(
    plot = signalPlot, label = format(peakPos, big.mark = ","), fontsize = 8,
    x = unit(peakPos, "native"), y = unit(1, "npc"),
    just = "bottom"
)

## ----arrange_example01, echo=FALSE, fig.height=4, fig.width=5, message=FALSE----
## Create page
bb_pageCreate(width = 3, height = 3, default.units = "inches")

## Draw arrows
bb_plotSegments(
    x0 = 0.5, y0 = 0, x1 = 0.5, y1 = 0.5,
    arrow = arrow(
        angle = 30,
        length = unit(0.1, "inches"),
        ends = "last",
        type = "closed"
    ),
    fill = "black"
)

bb_plotSegments(
    x0 = 0, y0 = 0.5, x1 = 0.5, y1 = 0.5,
    arrow = arrow(
        angle = 30,
        length = unit(0.1, "inches"),
        ends = "last",
        type = "closed"
    ),
    fill = "black"
)

## Draw text
bb_plotText(label = "0.5 in", x = 0.6, y = 0.25, just = "left")
bb_plotText(label = "0.5 in", x = 0.25, y = 0.6, just = "top")

## ----arrange_example02, echo=FALSE, fig.height=4, fig.width=5, message=FALSE----
## Create page
bb_pageCreate(width = 3, height = 3, default.units = "inches")

## Draw arrows
bb_plotSegments(
    x0 = 0.5, y0 = 0, x1 = 0.5, y1 = 0.5,
    arrow = arrow(
        angle = 30,
        length = unit(0.1, "inches"),
        ends = "last",
        type = "closed"
    ),
    fill = "black"
)

bb_plotSegments(
    x0 = 0, y0 = 0.5, x1 = 0.5, y1 = 0.5,
    arrow = arrow(
        angle = 30,
        length = unit(0.1, "inches"),
        ends = "last",
        type = "closed"
    ),
    fill = "black"
)

## Draw plot dimension arrows
bb_plotSegments(
    x0 = 0.5, y0 = 0.5, x1 = 2.5, y1 = 0.5,
    arrow = arrow(
        angle = 30,
        length = unit(0.1, "inches"),
        ends = "both",
        type = "closed"
    ),
    fill = "black"
)

bb_plotSegments(
    x0 = 0.5, y0 = 0.5, x1 = 0.5, y1 = 1.5,
    arrow = arrow(
        angle = 30,
        length = unit(0.1, "inches"),
        ends = "both",
        type = "closed"
    ),
    fill = "black"
)

## Draw text
bb_plotText(label = "2 in", x = 1.5, y = 0.4, just = "bottom")
bb_plotText(label = "1 in", x = 0.4, y = 1, just = "right")

## ----arrange_example03, echo=TRUE, fig.height=4, fig.width=5, message=FALSE----
## Create page
bb_pageCreate(width = 3, height = 3, default.units = "inches")

## Plot rectangle
bb_plotRect(
    x = 0.5, y = 0.5, width = 2, height = 1,
    just = c("left", "top"), default.units = "inches"
)

## ----arrange_example04, echo=TRUE, fig.height=4, fig.width=5, message=FALSE----
## Load data
data("bb_imrH3K27acData")

## Create page
bb_pageCreate(width = 3, height = 3, default.units = "inches")

## Define signal plot
signalPlot <- bb_plotSignal(
    data = bb_imrH3K27acData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    draw = FALSE
)

## Place plot on bb_page
bb_pagePlotPlace(
    plot = signalPlot,
    x = 0.5, y = 0.5, width = 2, height = 1,
    just = c("left", "top"), default.units = "inches"
)

## ----arrange_example05, echo=TRUE, fig.height=4, fig.width=5, message=FALSE----
# Load data
data("bb_imrH3K27acData")

## Create page
bb_pageCreate(width = 3, height = 3, default.units = "inches")

## Plot and place signal plot
signalPlot <- bb_plotSignal(
    data = bb_imrH3K27acData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 0.5, width = 2, height = 1,
    just = c("left", "top"), default.units = "inches"
)
## Remove signal plot
bb_pagePlotRemove(plot = signalPlot)

## ----arrange_example06, echo=FALSE, fig.height=5, fig.width=6, message=FALSE----
## Load libraries
library(BentoBox)

## Create page
bb_pageCreate(
    width = 4, height = 4, xgrid = 1, ygrid = 1,
    default.units = "inches"
)

## Define x and y coords, colors, and justifications
xcoords <- c(0.25, 2, 3.75, 3.75, 3.75, 2, 0.25, 0.25, 2)
ycoords <- c(0.25, 0.25, 0.25, 2, 3.75, 3.75, 3.75, 2, 2)
cols <- c(RColorBrewer::brewer.pal(n = 8, name = "Dark2"), "black")
justs <- list(
    c("left", "top"), c("center", "top"), c("right", "top"),
    c("right", "center"), c("right", "bottom"), c("center", "bottom"),
    c("left", "bottom"), c("left", "center"), c("center", "center")
)

## Define x and y shift vectors for adjusting text
s <- 0.1
xshift <- c(s, 0, -s, -s, -s, 0, s, s, 0)
yshift <- c(s, s, s, 0, -s, -s, -s, 0, -s * 3)

## Draw guides
bb_pageGuideHorizontal(y = ycoords, linecolor = "grey90", lty = 2)
bb_pageGuideVertical(x = xcoords, linecolor = "grey90", lty = 2)

## Draw box
bb_plotRect(
    x = 0.25, y = 0.25, width = 3.5, height = 3.5,
    just = c("left", "top")
)

## Plot points
bb_plotCircle(
    x = xcoords,
    y = ycoords,
    r = 0.035, linecolor = cols, fill = cols
)

## Text descriptions
for (i in 1:9) {
    bb_plotText(
        label = sprintf(
            "x = %s\ny = %s\njust = c('%s', '%s')",
            xcoords[i], ycoords[i], justs[[i]][1], justs[[i]][2]
        ),
        x = xcoords[i] + xshift[i],
        y = ycoords[i] + yshift[i],
        just = justs[[i]],
        fontsize = 6.5, fontcolor = cols[i], fontface = "bold"
    )
}

## ----arrange_example07, echo=FALSE, fig.height=5, fig.width=6, message=FALSE----
## Load libraries
library(BentoBox)

## Create page
bb_pageCreate(
    width = 4, height = 4, xgrid = 1, ygrid = 1,
    default.units = "inches"
)

## Define x and y coords, colors, and justifications
xcoords <- c(0.25, 2, 3.75, 3.75, 3.75, 2, 0.25, 0.25, 2)
ycoords <- c(0.25, 0.25, 0.25, 2, 3.75, 3.75, 3.75, 2, 2)
cols <- c(RColorBrewer::brewer.pal(n = 8, name = "Dark2"), "black")
justs <- list(
    c(0, 1), c(0.5, 1), c(1, 1),
    c(1, 0.5), c(1, 0), c(0.5, 0),
    c(0, 0), c(0, 0.5), c(0.5, 0.5)
)

## Define x and y shift vectors for adjusting text
s <- 0.1
xshift <- c(s, 0, -s, -s, -s, 0, s, s, 0)
yshift <- c(s, s, s, 0, -s, -s, -s, 0, -s * 3)

## Draw guides
bb_pageGuideHorizontal(y = ycoords, linecolor = "grey90", lty = 2)
bb_pageGuideVertical(x = xcoords, linecolor = "grey90", lty = 2)

## Draw box
bb_plotRect(
    x = 0.25, y = 0.25, width = 3.5, height = 3.5,
    just = c("left", "top")
)

## Plot points
bb_plotCircle(
    x = xcoords,
    y = ycoords,
    r = 0.035, linecolor = cols, fill = cols
)

## Text descriptions
for (i in 1:9) {
    bb_plotText(
        label = sprintf(
            "x = %s\ny = %s\njust = c(%s, %s)",
            xcoords[i], ycoords[i], justs[[i]][1], justs[[i]][2]
        ),
        x = xcoords[i] + xshift[i],
        y = ycoords[i] + yshift[i],
        just = justs[[i]],
        fontsize = 6.5, fontcolor = cols[i], fontface = "bold"
    )
}

## ----arrange_example08, echo=TRUE, fig.height=5, fig.width=6, message=FALSE----
## Load example Hi-C data
data("bb_imrHicData")

## Create a BentoBox page
bb_pageCreate(width = 3.25, height = 3.25, default.units = "inches")

## Plot Hi-C data with placing information
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
)

## Add color scale annotation with just = c("right", "top")
bb_annoHeatmapLegend(
    plot = hicPlot,
    x = 3, y = 0.25, width = 0.1, height = 1.25,
    just = c("right", "top"), default.units = "inches"
)

## ----plotting_example01, echo=TRUE, fig.height=6, fig.width=5, message=FALSE----
## Load example data
data("bb_imrHicData")
data("bb_bedpeData")
data("bb_imrH3K27acData")
data("bb_gwasData")

## Create a BentoBox page
bb_pageCreate(
    width = 3, height = 5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

## Plot Hi-C data in region
bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 0.5, width = 2, height = 2,
    just = c("left", "top"), default.units = "inches"
)

## Plot loop annotations
bb_plotPairsArches(
    data = bb_bedpeData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 2.5, width = 2, height = 0.25,
    just = c("left", "top"), default.units = "inches",
    fill = "black", linecolor = "black", flip = TRUE
)

## Plot signal track data
bb_plotSignal(
    data = bb_imrH3K27acData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 2.75, width = 2, height = 0.5,
    just = c("left", "top"), default.units = "inches"
)

## Plot GWAS data
bb_plotManhattan(
    data = bb_gwasData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    ymax = 1.1, cex = 0.20,
    x = 0.5, y = 3.5, width = 2, height = 0.5,
    just = c("left", "top"), default.units = "inches"
)

## Plot gene track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
bb_plotGenes(
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 4, width = 2, height = 0.5,
    just = c("left", "top"), default.units = "inches"
)

## Plot genome label
bb_plotGenomeLabel(
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 4.5, length = 2, scale = "Mb",
    just = c("left", "top"), default.units = "inches"
)

## ----eval=TRUE, message=FALSE, echo=TRUE--------------------------------------
params <- bb_params(
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, just = c("left", "top"),
    width = 2, length = 2, default.units = "inches"
)

## ----plotting_example02, echo=TRUE, fig.height=6, fig.width=5, message=FALSE----
## Load example data
data("bb_imrHicData")
data("bb_bedpeData")
data("bb_imrH3K27acData")
data("bb_gwasData")

## Create a BentoBox page
bb_pageCreate(
    width = 3, height = 5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

## Plot Hi-C data in region
bb_plotHicSquare(
    data = bb_imrHicData,
    params = params,
    y = 0.5, height = 2
)

## Plot loop annotations
bb_plotPairsArches(
    data = bb_bedpeData,
    params = params,
    y = 2.5, height = 0.25,
    fill = "black", linecolor = "black", flip = TRUE
)

## Plot signal track data
bb_plotSignal(
    data = bb_imrH3K27acData,
    params = params,
    y = 2.75, height = 0.5
)

## Plot GWAS data
bb_plotManhattan(
    data = bb_gwasData,
    params = params,
    ymax = 1.1, cex = 0.20,
    y = 3.5, height = 0.5
)

## Plot gene track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
bb_plotGenes(
    params = params,
    y = 4, height = 0.5
)

## Plot genome label
bb_plotGenomeLabel(
    params = params,
    y = 4.5, scale = "Mb"
)

## ----eval=TRUE, echo=TRUE, message=FALSE--------------------------------------
params <- bb_params(
    chrom = "chr21", chromstart = 29000000, chromend = 30000000,
    x = 0.5, just = c("left", "top"),
    width = 2, length = 2, default.units = "inches"
)

## ----plotting_example03, echo=FALSE, fig.height=6, fig.width=5, message=FALSE----
## Load example data
data("bb_imrHicData")
data("bb_bedpeData")
data("bb_imrH3K27acData")
data("bb_gwasData")

## Create a BentoBox page
bb_pageCreate(
    width = 3, height = 5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

## Plot Hi-C data in region
bb_plotHicSquare(
    data = bb_imrHicData,
    params = params,
    y = 0.5, height = 2
)

## Plot loop annotations
bb_plotPairsArches(
    data = bb_bedpeData,
    params = params,
    y = 2.5, height = 0.25,
    fill = "black", linecolor = "black", flip = TRUE
)

## Plot signal track data
bb_plotSignal(
    data = bb_imrH3K27acData,
    params = params,
    y = 2.75, height = 0.5
)

## Plot GWAS data
bb_plotManhattan(
    data = bb_gwasData,
    params = params,
    ymax = 1.1, cex = 0.20,
    y = 3.5, height = 0.5
)

## Plot gene track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
bb_plotGenes(
    params = params,
    y = 4, height = 0.5
)

## Plot genome label
bb_plotGenomeLabel(
    params = params,
    y = 4.5, scale = "Mb"
)

## ----eval=TRUE, echo=TRUE, message=FALSE--------------------------------------
params <- bb_params(
    gene = "LINC00113", geneBuffer = 100000, assembly = "hg19",
    x = 0.5, just = c("left", "top"),
    width = 2, length = 2, default.units = "inches"
)

## ----plotting_example04, echo=FALSE, fig.height=6, fig.width=5, message=FALSE----
## Load example data
data("bb_imrHicData")
data("bb_bedpeData")
data("bb_imrH3K27acData")
data("bb_gwasData")

## Create a BentoBox page
bb_pageCreate(
    width = 3, height = 5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

## Plot Hi-C data in region
bb_plotHicSquare(
    data = bb_imrHicData,
    params = params,
    y = 0.5, height = 2
)

## Plot loop annotations
bb_plotPairsArches(
    data = bb_bedpeData,
    params = params,
    y = 2.5, height = 0.25,
    fill = "black", linecolor = "black", flip = TRUE
)

## Plot signal track data
bb_plotSignal(
    data = bb_imrH3K27acData,
    params = params,
    y = 2.75, height = 0.5
)

## Plot GWAS data
bb_plotManhattan(
    data = bb_gwasData,
    params = params,
    ymax = 1.1, cex = 0.20,
    y = 3.5, height = 0.5
)

## Plot gene track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
bb_plotGenes(
    params = params,
    y = 4, height = 0.5
)

## Plot genome label
bb_plotGenomeLabel(
    params = params,
    y = 4.5, scale = "bp", fontsize = 7
)

## ----below_y, echo=TRUE, fig.height=6, fig.width=5, message=FALSE-------------
## Load example data
data("bb_imrHicData")
data("bb_bedpeData")
data("bb_imrH3K27acData")
data("bb_gwasData")

## bb_params
params <- bb_params(
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, just = c("left", "top"),
    width = 2, length = 2, default.units = "inches"
)

## Create a BentoBox page
bb_pageCreate(
    width = 3, height = 5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

## Plot Hi-C data in region
bb_plotHicSquare(
    data = bb_imrHicData,
    params = params,
    y = 0.5, height = 2
)

## Plot loop annotations
bb_plotPairsArches(
    data = bb_bedpeData,
    params = params,
    y = "0b",
    height = 0.25,
    fill = "black", linecolor = "black", flip = TRUE
)

## Plot signal track data
bb_plotSignal(
    data = bb_imrH3K27acData,
    params = params,
    y = "0b",
    height = 0.5
)

## Plot GWAS data
bb_plotManhattan(
    data = bb_gwasData,
    params = params,
    ymax = 1.1, cex = 0.20,
    y = "0.25b",
    height = 0.5
)

## Plot gene track
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
bb_plotGenes(
    params = params,
    y = "0b",
    height = 0.5
)

## Plot genome label
bb_plotGenomeLabel(
    params = params,
    y = "0b",
    scale = "Mb"
)

## ----ideogram_01, echo=TRUE, fig.height=5.25, fig.width=10.35, message=FALSE----
## Load cytoband data
data("cytoBand.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(GenomeInfoDb)

## Get sizes of chromosomes to scale their sizes
tx_db <- TxDb.Hsapiens.UCSC.hg19.knownGene
chromSizes <- GenomeInfoDb::seqlengths(tx_db)
maxChromSize <- max(chromSizes)

bb_pageCreate(
    width = 8.35, height = 4.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
xCoord <- 0.15
for (chr in c(paste0("chr", seq(1, 22)), "chrX", "chrY")) {
    height <- (4 * chromSizes[[chr]]) / maxChromSize
    bb_plotIdeogram(
        chrom = chr, assembly = "hg19",
        orientation = "v",
        x = xCoord, y = 4,
        width = 0.2, height = height,
        just = "bottom"
    )
    bb_plotText(
        label = gsub("chr", "", chr),
        x = xCoord, y = 4.1, fontsize = 10
    )
    xCoord <- xCoord + 0.35
}

## ----ideogram_02, echo=TRUE, fig.height=.5, fig.width=8.25, message=FALSE-----
bb_pageCreate(
    width = 6.25, height = 0.5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

bb_plotIdeogram(
    chrom = "chr1", assembly = "hg19",
    orientation = "h",
    x = 0.25, y = unit(0.25, "npc"), width = 5.75, height = 0.3
)

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  bb_plotIdeogram(
#      showBands = FALSE,
#      chrom = "chr1", assembly = "hg19",
#      orientation = "h",
#      x = 0.25, y = unit(0.25, "npc"), width = 5.75, height = 0.3
#  )

## ----ideogram_03, echo=FALSE, fig.height=.5, fig.width=8.25, message=FALSE----
bb_pageCreate(
    width = 6.25, height = 0.5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

bb_plotIdeogram(
    showBands = FALSE,
    chrom = "chr1", assembly = "hg19",
    orientation = "h",
    x = 0.25, y = unit(0.25, "npc"), width = 5.75, height = 0.3
)

## ----eval=TRUE, echo=FALSE, message=FALSE-------------------------------------
library(showtext)
font_add(
    family = "ProximaNova",
    regular = system.file("extdata",
        "proximanova-regular.otf",
        package = "BentoBox"
    ),
    bold = system.file("extdata",
        "proximanova-semibold.otf",
        package = "BentoBox"
    )
)
showtext_auto()

## ----ggplot_01, echo=TRUE, fig.showtext=TRUE, message=FALSE, warning=FALSE----
library(ggplot2)
library(scales)
data("bb_CasesUSA")

US_map <- ggplot(bb_CasesUSA, aes(long, lat, group = group)) +
    theme_void() +
    geom_polygon(aes(fill = cases_100K), color = "white", size = 0.3) +
    scale_fill_distiller(
        palette = "YlGnBu", direction = 1,
        labels = label_number(suffix = "", scale = 1e-3, accuracy = 1)
    ) +
    theme(
        legend.position = "left",
        legend.justification = c(0.5, 0.95),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.4, "cm"),
        plot.title = element_text(
            hjust = 0, vjust = -1,
            family = "ProximaNova", face = "bold",
            size = 12
        ),
        plot.title.position = "plot"
    ) +
    labs(title = "Thousands of COVID-19 Cases per 100,000 People") +
    coord_fixed(1.3)

print(US_map)

## ----ggplot_02, echo=TRUE, fig.showtext=TRUE, message=FALSE-------------------
data("bb_CasesNYFL")

# Format y-labels
ylabels <- seq(0, 2000000, by = 500000) / 1e6
ylabels[c(3, 5)] <- round(ylabels[c(3, 5)], digits = 0)
ylabels[c(2, 4)] <- round(ylabels[c(2, 4)], digits = 1)
ylabels[5] <- paste0(ylabels[5], "M cases")
ylabels[1] <- ""

bb_CasesNY <- bb_CasesNYFL[bb_CasesNYFL$state == "new york", ]
bb_CasesNYpoly <- rbind(
    bb_CasesNY,
    data.frame(
        "date" = as.Date("2021-03-07"),
        "state" = "new york",
        "caseIncrease" = -1 * sum(bb_CasesNY$caseIncrease)
    )
)

cases_NYline <- ggplot(
    bb_CasesNY,
    aes(x = date, y = cumsum(caseIncrease))
) +
    geom_polygon(data = bb_CasesNYpoly, fill = "#B8E6E6") +
    scale_x_date(
        labels = date_format("%b '%y"),
        breaks = as.Date(c("2020-05-01", "2020-09-01", "2021-01-01")),
        limits = as.Date(c("2020-01-29", "2021-03-07")),
        expand = c(0, 0)
    ) +
    scale_y_continuous(labels = ylabels, position = "right", expand = c(0, 0)) +
    geom_hline(
        yintercept = c(500000, 1000000, 1500000, 2000000),
        color = "white", linetype = "dashed", size = 0.3
    ) +
    coord_cartesian(ylim = c(0, 2000000)) +
    theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        text = element_text(family = "ProximaNova"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line.x.bottom = element_blank(),
        axis.line.y = element_line(size = 0.1, color = "#8F9BB3"),
        axis.text.x = element_text(
            size = 7, hjust = 0.5,
            vjust = 7.75, color = "black"
        ),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.2, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x.bottom = unit(-0.1, "cm"),
        plot.title = element_text(size = 8, hjust = 1),
        plot.title.position = "plot"
    )

print(cases_NYline)

## ----ggplot_03, echo=TRUE, fig.showtext=TRUE, message=FALSE-------------------
bb_CasesFL <- bb_CasesNYFL[bb_CasesNYFL$state == "florida", ]
bb_CasesFLpoly <- rbind(
    bb_CasesFL,
    data.frame(
        "date" = as.Date("2021-03-07"),
        "state" = "florida",
        "caseIncrease" = -1 * sum(bb_CasesFL$caseIncrease)
    )
)

cases_FLline <- ggplot(
    bb_CasesFL,
    aes(x = date, y = cumsum(caseIncrease))
) +
    geom_polygon(data = bb_CasesFLpoly, fill = "#B8E6E6") +
    scale_x_date(
        labels = date_format("%b '%y"),
        breaks = as.Date(c("2020-05-01", "2020-09-01", "2021-01-01")),
        limits = as.Date(c("2020-01-29", "2021-03-07")),
        expand = c(0, 0)
    ) +
    scale_y_continuous(labels = ylabels, position = "right", expand = c(0, 0)) +
    geom_hline(
        yintercept = c(500000, 1000000, 1500000, 2000000),
        color = "white", linetype = "dashed", size = 0.3
    ) +
    coord_cartesian(ylim = c(0, 2000000)) +
    theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        text = element_text(family = "ProximaNova"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.line.y = element_line(size = 0.1, color = "#8F9BB3"),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.text.x = element_text(
            size = 7, hjust = 0.5,
            vjust = 7.75, color = "black"
        ),
        axis.ticks = element_line(color = "black", size = 0.2),
        axis.ticks.y = element_blank(),
        axis.ticks.length.x.bottom = unit(-0.1, "cm"),
        plot.title = element_text(size = 8, hjust = 1),
        plot.title.position = "plot"
    )

print(cases_FLline)

## ----ggplot_04, echo=TRUE, fig.showtext=TRUE, message=FALSE-------------------
data("bb_VaccinesNYFL")

vaccines_NYpie <- ggplot(
    bb_VaccinesNYFL[bb_VaccinesNYFL$state == "new york", ],
    aes(x = "", y = value, fill = vax_group)
) +
    geom_bar(width = 1, stat = "identity") +
    theme_void() +
    scale_fill_manual(values = c("#FBAA7E", "#F7EEBF", "#FBCB88")) +
    coord_polar(theta = "y", start = 2.125, clip = "off") +
    geom_text(aes(
        x = c(1.9, 2, 1.9),
        y = c(1.65e7, 1.3e6, 7.8e6),
        label = paste0(percent, "%")
    ),
    size = 2.25, color = "black"
    ) +
    theme(
        legend.position = "none",
        plot.title = element_text(
            hjust = 0.5, vjust = -3.5, size = 10,
            family = "ProximaNova", face = "bold"
        ),
        text = element_text(family = "ProximaNova")
    ) +
    labs(title = "New York")

print(vaccines_NYpie)

## ----ggplot_05, echo = TRUE, eval = TRUE--------------------------------------
vaccines_FLpie <- ggplot(
    bb_VaccinesNYFL[bb_VaccinesNYFL$state == "florida", ],
    aes(x = "", y = value, fill = vax_group)
) +
    geom_bar(width = 1, stat = "identity") +
    scale_fill_manual(values = c("#FBAA7E", "#F7EEBF", "#FBCB88")) +
    theme_void() +
    coord_polar(theta = "y", start = pi / 1.78, clip = "off") +
    geom_text(aes(
        x = c(1.95, 2, 1.9),
        y = c(1.9e7, 1.83e6, 9.6e6),
        label = paste0(percent, "%")
    ),
    color = "black",
    size = 2.25
    ) +
    theme(
        legend.position = "none",
        plot.title = element_text(
            hjust = 0.5, vjust = -4, size = 10,
            family = "ProximaNova", face = "bold"
        ),
        text = element_text(family = "ProximaNova")
    ) +
    labs(title = "Florida")

print(vaccines_FLpie)

## ----gg_plot05, echo=TRUE, fig.height=4.5, fig.showtext=TRUE, fig.width=11.5----
bb_pageCreate(width = 9.5, height = 3.5, default.units = "inches")

bb_plotGG(
    plot = US_map,
    x = 0.1, y = 0,
    width = 6.5, height = 3.5, just = c("left", "top")
)
bb_plotGG(
    plot = cases_NYline,
    x = 6.25, y = 1.8,
    width = 3.025, height = 1.4, just = c("left", "bottom")
)
bb_plotGG(
    plot = cases_FLline,
    x = 6.25, y = 3.5,
    width = 3.025, height = 1.4, just = c("left", "bottom")
)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  bb_plotGG(
#      plot = vaccines_NYpie,
#      x = 6.37, y = -0.05,
#      width = 1.45, height = 1.45, just = c("left", "top")
#  )
#  bb_plotGG(
#      plot = vaccines_FLpie,
#      x = 6.37, y = 1.67,
#      width = 1.45, height = 1.45, just = c("left", "top")
#  )

## ----gg_plot06, echo=FALSE, fig.height=4.5, fig.showtext=TRUE, ,fig.width=11.5----
bb_pageCreate(width = 9.5, height = 3.5, default.units = "inches")

bb_plotGG(
    plot = US_map,
    x = 0.1, y = 0,
    width = 6.5, height = 3.5, just = c("left", "top")
)
bb_plotGG(
    plot = cases_NYline,
    x = 6.25, y = 1.8,
    width = 3.025, height = 1.4, just = c("left", "bottom")
)
bb_plotGG(
    plot = cases_FLline,
    x = 6.25, y = 3.5,
    width = 3.025, height = 1.4, just = c("left", "bottom")
)
bb_plotGG(
    plot = vaccines_NYpie,
    x = 6.37, y = -0.05,
    width = 1.45, height = 1.45, just = c("left", "top")
)
bb_plotGG(
    plot = vaccines_FLpie,
    x = 6.37, y = 1.67,
    width = 1.45, height = 1.45, just = c("left", "top")
)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  bb_plotText(
#      label = c("not", "partially", "fully vaccinated"),
#      fontfamily = "ProximaNova", fontcolor = "black", fontsize = 7,
#      x = c(6.58, 7.3, 7.435),
#      y = c(0.74, 1.12, 0.51), just = c("left", "bottom")
#  )
#  bb_plotText(
#      label = c("not", "partially", "fully vaccinated"),
#      fontfamily = "ProximaNova", fontcolor = "black", fontsize = 7,
#      x = c(6.58, 7.39, 7.435),
#      y = c(2.47, 2.75, 2.2), just = c("left", "bottom")
#  )

## ----gg_plot07, echo=FALSE, fig.height=4.5, fig.showtext=TRUE, fig.width=11.5----
bb_pageCreate(
    width = 9.5, height = 3.5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
bb_plotGG(
    plot = US_map,
    x = 0.1, y = 0,
    width = 6.5, height = 3.5, just = c("left", "top")
)
bb_plotGG(
    plot = cases_NYline,
    x = 6.25, y = 1.8,
    width = 3.025, height = 1.4, just = c("left", "bottom")
)
bb_plotGG(
    plot = cases_FLline,
    x = 6.25, y = 3.5,
    width = 3.025, height = 1.4, just = c("left", "bottom")
)
bb_plotGG(
    plot = vaccines_NYpie,
    x = 6.37, y = -0.05,
    width = 1.45, height = 1.45, just = c("left", "top")
)
bb_plotGG(
    plot = vaccines_FLpie,
    x = 6.37, y = 1.67,
    width = 1.45, height = 1.45, just = c("left", "top")
)
bb_plotText(
    label = c("not", "partially", "fully vaccinated"),
    fontfamily = "ProximaNova", fontcolor = "black", fontsize = 7,
    x = c(6.58, 7.3, 7.435),
    y = c(0.74, 1.12, 0.51), just = c("left", "bottom")
)
bb_plotText(
    label = c("not", "partially", "fully vaccinated"),
    fontfamily = "ProximaNova", fontcolor = "black", fontsize = 7,
    x = c(6.58, 7.39, 7.435),
    y = c(2.47, 2.75, 2.2), just = c("left", "bottom")
)

## ----shapes, echo=TRUE, fig.height=7, fig.showtext=TRUE, fig.width=6----------
library(png)
library(showtext)
font_add(
    family = "ProximaNova",
    regular = system.file("extdata",
        "proximanova-regular.otf",
        package = "BentoBox"
    ),
    bold = system.file("extdata",
        "proximanova-semibold.otf",
        package = "BentoBox"
    )
)
showtext_auto()

edamaman <- readPNG(system.file("images",
    "bento-edamaman.png",
    package = "BentoBox"
))
logotype <- readPNG(system.file("images",
    "bento-logotype-singleline-black.png",
    package = "BentoBox"
))

bb_pageCreate(
    width = 5, height = 6, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
bb_plotRaster(
    image = logotype,
    x = 2.5, y = 0.25, width = 3.25, height = 0.5, just = "top"
)
bb_plotRaster(
    image = edamaman,
    x = 2.5, y = 5.5, width = 2, height = 4, just = "bottom"
)
bb_plotText(
    label = "Edamaman",
    fontsize = 20, fontfamily = "ProximaNova", fontface = "bold",
    x = 2.5, y = 0.9, just = "top"
)

## ----anno_01, echo=TRUE, eval=TRUE, fig.height=5, fig.width=6-----------------
data("bb_imrHicData")
bb_pageCreate(
    width = 3, height = 3.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
)

bb_annoGenomeLabel(
    plot = hicPlot, scale = "Mb",
    x = 0.25, y = 2.76
)

## ----anno_02, echo=FALSE, eval=TRUE, warning=TRUE, fig.height=5, fig.width=6----
data("bb_imrHicData")
bb_pageCreate(
    width = 3, height = 3.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28255554, chromend = 29354665,
    x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
)
bb_annoGenomeLabel(
    plot = hicPlot, scale = "Mb",
    x = 0.25, y = 2.76
)

## ----anno_03, echo=TRUE, eval=TRUE, message=FALSE, fig.height=5, fig.width=6----
data("bb_imrHicData")
bb_pageCreate(
    width = 3, height = 3.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28255554, chromend = 29354665,
    x = 0.25, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
)
bb_annoGenomeLabel(
    plot = hicPlot, scale = "bp",
    x = 0.25, y = 2.76
)

## ----anno_04, echo=FALSE, fig.height=2.5, fig.width=7, message=FALSE----------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg19)

bb_pageCreate(
    width = 5, height = 1.5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
geneTrack <- bb_plotGenes(
    chrom = "chr1", chromstart = 110146000, chromend = 110146050,
    x = 0.5, y = 0.25, width = 4, height = 0.75, default.units = "inches"
)
bb_annoGenomeLabel(plot = geneTrack, x = 0.5, y = 1.01)

## ----anno_05, echo=FALSE, fig.height=2.5, fig.width=7, message=FALSE----------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

bb_pageCreate(
    width = 5, height = 1.5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
geneTrack <- bb_plotGenes(
    chrom = "chr1", chromstart = 110146000, chromend = 110146025,
    x = 0.5, y = 0.25, width = 4, height = 0.75, default.units = "inches"
)
bb_annoGenomeLabel(plot = geneTrack, x = 0.5, y = 1.01)

## ----anno_06, echo=TRUE, eval=TRUE, message=FALSE, fig.height=5, fig.width=6----
data("bb_imrHicData")
bb_pageCreate(
    width = 3.25, height = 3, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 0.25, width = 2.5, height = 2.5, default.units = "inches"
)
bb_annoGenomeLabel(
    plot = hicPlot, scale = "Mb",
    axis = "y",
    x = 0.5, y = 0.25,
    just = c("right", "top")
)

## ----anno_07, echo=TRUE, fig.height=3.75, fig.width=9.5, message=FALSE--------
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
data("bb_gwasData")

bb_pageCreate(
    width = 7.5, height = 2.75, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
manhattanPlot <- bb_plotManhattan(
    data = bb_gwasData, assembly = "hg19",
    fill = c("grey", "#37a7db"),
    sigLine = TRUE,
    col = "grey", lty = 2, range = c(0, 14),
    x = 0.5, y = 0.25, width = 6.5, height = 2,
    just = c("left", "top"),
    default.units = "inches"
)
bb_annoGenomeLabel(
    plot = manhattanPlot, x = 0.5, y = 2.25, fontsize = 8,
    just = c("left", "top"), default.units = "inches"
)
bb_plotText(
    label = "Chromosome", fontsize = 8,
    x = 3.75, y = 2.45, just = "center", default.units = "inches"
)

## Annotate y-axis
bb_annoYaxis(
    plot = manhattanPlot, at = c(0, 2, 4, 6, 8, 10, 12, 14),
    axisLine = TRUE, fontsize = 8
)
## Plot y-axis label
bb_plotText(
    label = "-log10(p-value)", x = 0.15, y = 1.25, rot = 90,
    fontsize = 8, fontface = "bold", just = "center",
    default.units = "inches"
)

## ----anno_08, echo=TRUE, fig.height=4.25, fig.width=5.25, message=FALSE-------
data("bb_imrHicData")

bb_pageCreate(
    width = 3.25, height = 3.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
params <- bb_params(
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    assembly = "hg19",
    x = 0.25, width = 2.75, just = c("left", "top"), default.units = "inches"
)
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData, params = params,
    zrange = c(0, 70), resolution = 10000,
    y = 0.25, height = 2.75
)

## Annotate Hi-C heatmap legend
bb_annoHeatmapLegend(
    plot = hicPlot, fontsize = 7,
    orientation = "v",
    x = 0.125, y = 0.25,
    width = 0.07, height = 0.5, just = c("left", "top"),
    default.units = "inches"
)

bb_annoHeatmapLegend(
    plot = hicPlot, fontsize = 7,
    orientation = "h",
    x = 3, y = 3.055,
    width = 0.5, height = 0.07, just = c("right", "top"),
    default.units = "inches"
)

## ----anno_09, echo=TRUE, fig.height=4.25, fig.width=5.25, message=FALSE-------
data("bb_imrHicData")
data("bb_bedpeData")

bb_pageCreate(
    width = 3.25, height = 3.24, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData, resolution = 10000, zrange = c(0, 70),
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.25, y = 0.25, width = 2.75, height = 2.75,
    just = c("left", "top"),
    default.units = "inches"
)

## Annotate pixels
pixels <- bb_annoPixels(
    plot = hicPlot, data = bb_bedpeData, type = "box",
    half = "top"
)

## ----anno_10, echo=TRUE, fig.height=4.25, fig.width=5.25, message=FALSE-------
data("bb_imrHicData")
data("bb_bedpeData")

## Subset BEDPE data
bb_bedpeData <- bb_bedpeData[which(bb_bedpeData$start1 == 28220000 &
    bb_bedpeData$start2 == 29070000), ]

bb_pageCreate(
    width = 3.25, height = 3.24, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData, resolution = 10000, zrange = c(0, 70),
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.25, y = 0.25, width = 2.75, height = 2.75,
    just = c("left", "top"),
    default.units = "inches"
)

## Annotate pixel
pixels <- bb_annoPixels(
    plot = hicPlot, data = bb_bedpeData, type = "arrow",
    half = "bottom", shift = 12
)

## ----anno_11, echo=TRUE, fig.height=3.25, fig.width=8.25, message=FALSE-------
data("cytoBand.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

bb_pageCreate(
    width = 6.25, height = 2.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
ideogramPlot <- bb_plotIdeogram(
    chrom = "chr21", assembly = "hg19",
    orientation = "h",
    x = 0.25, y = 0.5, width = 5.75, height = 0.3, just = "left"
)

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  region <- bb_params(chrom = "chr21", chromstart = 28000000, chromend = 30300000)
#  bb_annoHighlight(
#      plot = ideogramPlot, params = region,
#      fill = "red",
#      y = 0.25, height = 0.5, just = c("left", "top"), default.units = "inches"
#  )

## ----anno_12, echo=FALSE, fig.height=3.25, fig.width=8.25, message=FALSE------
data("cytoBand.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
bb_pageCreate(
    width = 6.25, height = 2.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

ideogramPlot <- bb_plotIdeogram(
    chrom = "chr21", assembly = "hg19",
    orientation = "h",
    x = 0.25, y = 0.5, width = 5.75, height = 0.3, just = "left"
)
region <- bb_params(chrom = "chr21", chromstart = 28000000, chromend = 30300000)

bb_annoHighlight(
    plot = ideogramPlot, params = region,
    fill = "red",
    y = 0.25, height = 0.5, just = c("left", "top"), default.units = "inches"
)

## ----eval=FALSE, echo=TRUE, message=FALSE-------------------------------------
#  bb_annoZoomLines(
#      plot = ideogramPlot, params = region,
#      y0 = 0.75, x1 = c(0.25, 6), y1 = 1.25, default.units = "inches"
#  )

## ----anno_13, echo=FALSE, fig.height=3.25, fig.width=8.25, message=FALSE------
data("cytoBand.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
bb_pageCreate(
    width = 6.25, height = 2.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

ideogramPlot <- bb_plotIdeogram(
    chrom = "chr21", assembly = "hg19",
    orientation = "h",
    x = 0.25, y = 0.5, width = 5.75, height = 0.3, just = "left"
)
region <- bb_params(chrom = "chr21", chromstart = 28000000, chromend = 30300000)

bb_annoHighlight(
    plot = ideogramPlot, params = region,
    fill = "red",
    y = 0.25, height = 0.5, just = c("left", "top"), default.units = "inches"
)
bb_annoZoomLines(
    plot = ideogramPlot, params = region,
    y0 = 0.75, x1 = c(0.25, 6), y1 = 1.25, default.units = "inches"
)

## ----anno_14, echo=FALSE, fig.height=3.25, fig.width=8.25, message=FALSE------
data("cytoBand.Hsapiens.UCSC.hg19")
data("bb_imrH3K27acData")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
bb_pageCreate(
    width = 6.25, height = 2.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

ideogramPlot <- bb_plotIdeogram(
    chrom = "chr21", assembly = "hg19",
    orientation = "h",
    x = 0.25, y = 0.5, width = 5.75, height = 0.3, just = "left"
)
region <- bb_params(chrom = "chr21", chromstart = 28000000, chromend = 30300000)

bb_annoHighlight(
    plot = ideogramPlot, params = region,
    fill = "red",
    y = 0.25, height = 0.5, just = c("left", "top"), default.units = "inches"
)
bb_annoZoomLines(
    plot = ideogramPlot, params = region,
    y0 = 0.75, x1 = c(0.25, 6), y1 = 1.25, default.units = "inches"
)


signalPlot <- bb_plotSignal(
    data = bb_imrH3K27acData, params = region,
    x = 0.25, y = 1.25, width = 5.75, height = 0.65
)

bb_annoGenomeLabel(
    plot = signalPlot, scale = "Mb",
    x = 0.25, y = 1.95
)

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
bb_defaultPackages("hg19")
bb_defaultPackages("hg38")
bb_defaultPackages("mm10")

## ----echo=TRUE, eval=TRUE-----------------------------------------------------
bb_genomes()

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  library(GenomicFeatures)
#  TxDb.Hsapiens.Ensembl.GrCh38.103 <- makeTxDbFromEnsembl(
#      organism =
#          "Homo sapiens"
#  )

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  Ensembl38 <- bb_assembly(
#      Genome = "Ensembl.GRCh38.103",
#      TxDb = TxDb.Hsapiens.Ensembl.GrCh38.103,
#      OrgDb = "org.Hs.eg.db",
#      BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38",
#      gene.id = "ENSEMBL", display.column = "SYMBOL"
#  )

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  bb_plotGenes(
#      chrom = "chr8", chromstart = 1000000, chromend = 2000000,
#      bg = "#f6f6f6",
#      x = 0.5, y = 0.5, width = 2, height = 1, just = c("left", "top"),
#      default.units = "inches"
#  )

## ----aes_01, eval=TRUE, echo=FALSE, message=FALSE, fig.width=5, fig.height=3----
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
bb_pageCreate(
    width = 3, height = 2, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
genomicRegion <- bb_params(
    chrom = "chr8",
    chromstart = 1000000, chromend = 2000000
)
genesPlot <- bb_plotGenes(
    params = genomicRegion, bg = "#f6f6f6",
    x = 0.5, y = 0.5, width = 2, height = 1,
    just = c("left", "top"), default.units = "inches"
)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  bb_plotRanges(
#      data = bb_bedData,
#      chrom = "chr21", chromstart = 29073000, chromend = 29074000,
#      fill = c("#7ecdbb", "#37a7db"),
#      baseline = TRUE, baseline.color = "black",
#      x = 0.5, y = 0.25, width = 6.5, height = 4.25,
#      just = c("left", "top"), default.units = "inches"
#  )

## ----aes_02, echo=FALSE, fig.height=6, fig.width=9.5--------------------------
data("bb_bedData")
bb_pageCreate(
    width = 7.5, height = 5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

pileupPlot <- bb_plotRanges(
    data = bb_bedData,
    chrom = "chr21", chromstart = 29073000, chromend = 29074000,
    fill = c("#7ecdbb", "#37a7db"),
    x = 0.5, y = 0.25, width = 6.5, height = 4.25,
    baseline = TRUE, baseline.color = "black",
    just = c("left", "top"), default.units = "inches"
)

## ----eval = TRUE, echo = TRUE-------------------------------------------------
data("bb_bedData")
head(bb_bedData)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  bb_plotRanges(
#      data = bb_bedData,
#      chrom = "chr21", chromstart = 29073000, chromend = 29074000,
#      fill = c("#7ecdbb", "#37a7db"),
#      colorby = colorby("strand"),
#      x = 0.5, y = 0.25, width = 6.5, height = 4.25,
#      just = c("left", "top"), default.units = "inches"
#  )

## ----aes_03, echo=FALSE, fig.height=6, fig.width=9.5--------------------------

bb_pageCreate(
    width = 7.5, height = 5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)

bb_plotRanges(
    data = bb_bedData,
    chrom = "chr21", chromstart = 29073000, chromend = 29074000,
    fill = c("#7ecdbb", "#37a7db"),
    colorby = colorby("strand"),
    strandSplit = TRUE,
    x = 0.5, y = 0.25, width = 6.5, height = 4.25,
    just = c("left", "top"), default.units = "inches"
)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
data("bb_bedData")
bb_bedData$strand <- as.factor(bb_bedData$strand)
levels(bb_bedData$strand) <- c("+", "-")
head(bb_bedData$strand)

## ----aes_04, echo=FALSE, fig.height=6, fig.width=9.5--------------------------

bb_pageCreate(
    width = 7.5, height = 5, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
bb_plotRanges(
    data = bb_bedData,
    chrom = "chr21", chromstart = 29073000, chromend = 29074000,
    fill = c("#7ecdbb", "#37a7db"),
    colorby = colorby("strand"),
    strandSplit = TRUE,
    x = 0.5, y = 0.25, width = 6.5, height = 4.25,
    just = c("left", "top"), default.units = "inches"
)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
data("bb_bedpeData")
bb_bedpeData$length <- (bb_bedpeData$start2 - bb_bedpeData$start1) / 1000
head(bb_bedpeData$length)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  bedpePlot <- bb_plotPairsArches(
#      data = bb_bedpeData,
#      chrom = "chr21", chromstart = 27900000, chromend = 30700000,
#      fill = colorRampPalette(c("dodgerblue2", "firebrick2")),
#      linecolor = "fill",
#      colorby = colorby("length"),
#      archHeight = bb_bedpeData$length / max(bb_bedpeData$length),
#      alpha = 1,
#      x = 0.25, y = 0.25, width = 7, height = 1.5,
#      just = c("left", "top"),
#      default.units = "inches"
#  )

## ----aes_05, echo=FALSE, fig.width=9.5, message=FALSE, height=3---------------

bb_pageCreate(
    width = 7.5, height = 2, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
heights <- bb_bedpeData$length / max(bb_bedpeData$length)
bedpePlot <- bb_plotPairsArches(
    data = bb_bedpeData,
    chrom = "chr21", chromstart = 27900000, chromend = 30700000,
    fill = colorRampPalette(c("dodgerblue2", "firebrick2")),
    linecolor = "fill",
    colorby = colorby("length"), archHeight = heights, alpha = 1,
    x = 0.25, y = 0.25, width = 7, height = 1.5,
    just = c("left", "top"),
    default.units = "inches"
)

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  bb_annoHeatmapLegend(
#      plot = bedpePlot, fontcolor = "black",
#      x = 7.0, y = 0.25,
#      width = 0.10, height = 1, fontsize = 10
#  )

## ----aes_06, echo=FALSE, fig.width=9.5, message=FALSE, height=3---------------

bb_pageCreate(
    width = 7.5, height = 2, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
heights <- bb_bedpeData$length / max(bb_bedpeData$length)
bedpePlot <- bb_plotPairsArches(
    data = bb_bedpeData,
    chrom = "chr21", chromstart = 27900000, chromend = 30700000,
    fill = colorRampPalette(c("dodgerblue2", "firebrick2")),
    linecolor = "fill",
    colorby = colorby("length"), archHeight = heights, alpha = 1,
    x = 0.25, y = 0.25, width = 7, height = 1.5,
    just = c("left", "top"),
    default.units = "inches"
)
bb_annoHeatmapLegend(
    plot = bedpePlot, fontcolor = "black",
    x = 7.0, y = 0.25,
    width = 0.10, height = 1, fontsize = 10
)

## ----aes_07, eval=TRUE, echo=TRUE, message=FALSE, fig.width=5, fig.height=1.5----
bb_pageCreate(
    width = 5, height = 1.25,
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
genePlot <- bb_plotGenes(
    chrom = "chr2", chromstart = 1000000, chromend = 20000000,
    x = 0.25, y = 0.25, width = 4.75, height = 1
)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
genePlot$genes

## ----aes_08, echo=TRUE, fig.height=1.5, fig.width=5, message=FALSE------------
bb_pageCreate(
    width = 5, height = 1.25,
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
genePlot <- bb_plotGenes(
    chrom = "chr2", chromstart = 1000000, chromend = 20000000,
    geneOrder = c("MIR3125"),
    x = 0.25, y = 0.25, width = 4.75, height = 1
)

## ----aes_09, echo=TRUE, fig.height=2.25, fig.width=5, message=FALSE-----------
new_hg19 <- bb_assembly(
    Genome = "id_hg19",
    TxDb = "TxDb.Hsapiens.UCSC.hg19.knownGene",
    OrgDb = "org.Hs.eg.db",
    gene.id.column = "ENTREZID",
    display.column = "ENSEMBL"
)
bb_pageCreate(
    width = 5, height = 1.25,
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
genePlot <- bb_plotGenes(
    chrom = "chr2", chromstart = 1000000, chromend = 20000000,
    assembly = new_hg19,
    x = 0.25, y = 0.25, width = 4.75, height = 1
)

## ----aes_10, echo=TRUE, fig.height=6, fig.width=7-----------------------------
bb_pageCreate(
    width = 6, height = 5,
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
transcriptPlot <- bb_plotTranscripts(
    chrom = "chr2", chromstart = 1000000, chromend = 20000000,
    labels = "both",
    x = 0.25, y = 0.25, width = 5.5, height = 4.5
)

## ----eval=TRUE, echo=TRUE, message=FALSE--------------------------------------
geneHighlights <- data.frame("geneName" = "RRM2", "color" = "steel blue")

## ----aes_11, eval=TRUE, echo=TRUE, message=FALSE, fig.width=5, fig.height=1.5----
bb_pageCreate(
    width = 5, height = 1.25,
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
genePlot <- bb_plotGenes(
    chrom = "chr2", chromstart = 1000000, chromend = 20000000,
    geneHighlights = geneHighlights, geneBackground = "grey",
    x = 0.25, y = 0.25, width = 4.75, height = 1
)

## ----aes_12, eval=TRUE, echo=TRUE, message=FALSE, fig.width=5, fig.height=1.5----
geneHighlights <- data.frame(
    "geneName" = c("RRM2", "PXDN"),
    "color" = c("steel blue", "red")
)

bb_pageCreate(
    width = 5, height = 1.25,
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
genePlot <- bb_plotGenes(
    chrom = "chr2", chromstart = 1000000, chromend = 20000000,
    geneHighlights = geneHighlights, geneBackground = "grey",
    x = 0.25, y = 0.25, width = 4.75, height = 1
)

## ----aes_13, echo=TRUE, fig.height=6, fig.width=7-----------------------------
bb_pageCreate(
    width = 6, height = 5,
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
transcriptPlot <- bb_plotTranscripts(
    chrom = "chr2", chromstart = 1000000, chromend = 20000000,
    strandSplit = TRUE,
    x = 0.25, y = 0.25, width = 5.5, height = 4.5
)

## ----aes_14, echo=FALSE, fig.height=3.9, fig.width=9.75, message=FALSE--------
data("bb_imrHicData")

bb_pageCreate(
    width = 7.1, height = 2.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
params <- bb_params(
    chrom = "chr21", chromstart = 28950000, chromend = 29800000,
    assembly = "hg19", resolution = 10000
)

hicPlot_square <- bb_plotHicSquare(
    data = bb_imrHicData, params = params,
    zrange = c(0, 70),
    x = 0.25, y = 0.25,
    width = 1.6, height = 1.6
)
bb_annoGenomeLabel(
    plot = hicPlot_square, fontsize = 8,
    x = 0.25, y = 1.9, scale = "Kb"
)
bb_annoGenomeLabel(
    plot = hicPlot_square, axis = "y", fontsize = 8,
    x = 0.20, y = 0.25, scale = "Kb",
    just = c("right", "top")
)
hicPlot_triangle <- bb_plotHicTriangle(
    data = bb_imrHicData, params = params,
    zrange = c(0, 70),
    x = 2.1, y = 1.85,
    width = 2.25, height = 2,
    just = c("left", "bottom")
)
bb_annoGenomeLabel(
    plot = hicPlot_triangle, fontsize = 8,
    x = 2.1, y = 1.9, scale = "Kb"
)
hicPlot_rectangle <- bb_plotHicRectangle(
    data = bb_imrHicData, params = params,
    zrange = c(0, 70),
    x = 4.6, y = 1.85,
    width = 2.25, height = 1.125,
    just = c("left", "bottom")
)
bb_annoGenomeLabel(
    plot = hicPlot_rectangle, fontsize = 8,
    x = 4.6, y = 1.9, scale = "Kb"
)

## ----aes_15, echo=TRUE, fig.height=4.25, fig.width=5.25, message=FALSE--------
data("bb_gmHicData")
data("bb_imrHicData")

bb_pageCreate(
    width = 3.25, height = 3.25, default.units = "inches",
    showGuides = FALSE, xgrid = 0, ygrid = 0
)
params <- bb_params(
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    assembly = "hg19", resolution = 10000,
    x = 0.25, width = 2.75, just = c("left", "top"), default.units = "inches"
)


hicPlot_top <- bb_plotHicSquare(
    data = bb_gmHicData, params = params,
    zrange = c(0, 200),
    half = "top",
    y = 0.25, height = 2.75
)
hicPlot_bottom <- bb_plotHicSquare(
    data = bb_imrHicData, params = params,
    zrange = c(0, 70),
    half = "bottom",
    y = 0.25, height = 2.75
)

