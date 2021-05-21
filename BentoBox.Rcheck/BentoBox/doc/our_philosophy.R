## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.height = 5,
    fig.width = 6,
    fig.align = "center"
)
library(BentoBox)
library(grid)

## ----philosophy_1, eval=TRUE--------------------------------------------------
bb_pageCreate(width = 3, height = 3, default.units = "inches")

## ----philosphy_2, echo=FALSE, eval=TRUE, message=FALSE------------------------
bb_pageCreate(width = 3, height = 3, default.units = "inches")
bb_plotSegments(x0 = 0.5, y0 = 0, x1 = 0.5, y1 = 0.5, 
                arrow = arrow(angle = 30, 
                            length = unit(0.1, "inches"), 
                            ends = "last", type = "closed"), 
                fill = "black")
bb_plotSegments(x0 = 0, y0 = 0.5, x1 = 0.5, y1 = 0.5, 
                arrow = arrow(angle = 30, 
                            length = unit(0.1, "inches"), 
                            ends = "last", type = "closed"),
                fill = "black")
bb_plotText(label = "0.5 in", x = 0.6, y = 0.25, just = "left")
bb_plotText(label = "0.5 in", x = 0.25, y = 0.6, just = "top")

## ----philosphy_3, echo=FALSE, eval=TRUE, message=FALSE------------------------
bb_pageCreate(width = 3, height = 3, default.units = "inches")
bb_plotSegments(x0 = 0.5, y0 = 0, x1 = 0.5, y1 = 0.5, 
                arrow = arrow(angle = 30, 
                            length = unit(0.1, "inches"), 
                            ends = "last", type = "closed"), 
                fill = "black")
bb_plotSegments(x0 = 0, y0 = 0.5, x1 = 0.5, y1 = 0.5, 
                arrow = arrow(angle = 30, 
                            length = unit(0.1, "inches"), 
                            ends = "last", type = "closed"), 
                fill = "black")
bb_plotSegments(x0 = 0.5, y0 = 0.5, x1 = 2.5, y1 = 0.5, fill = "black", 
                arrow = arrow(angle = 30, 
                            length = unit(0.1, "inches"), 
                            ends = "both", type = "closed"))
bb_plotSegments(x0 = 0.5, y0 = 0.5, x1 = 0.5, y1 = 1.5, fill = "black", 
                arrow = arrow(angle = 30, 
                            length = unit(0.1, "inches"), 
                            ends = "both", type = "closed"))
bb_plotText(label = "2 in", x = 1.5, y = 0.4, just = "bottom")
bb_plotText(label = "1 in", x = 0.4, y = 1, just = "right")

## ----philosophy_4, echo=TRUE, eval=TRUE, message=FALSE------------------------
## Create BentoBox page
bb_pageCreate(width = 3, height = 3, default.units = "inches")

## Load signal data
data("bb_imrH3K27acData")

## Plot and place signal data with precise measurements
signalPlot <- bb_plotSignal(
    data = bb_imrH3K27acData,
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = 0.5, y = 0.5, width = 2, height = 1,
    just = c("left", "top"), default.units = "inches"
)

## ----philosphy_5, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE----------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

## Create BentoBox page
bb_pageCreate(width = 3, height = 2, default.units = "inches")

## Define genomic region with `bb_params`
genomicRegion <- bb_params(chrom = "chr8", 
                        chromstart = 1000000, chromend = 2000000)

## Plot genes with background color highlighting size
genesPlot <- bb_plotGenes(
    params = genomicRegion, bg = "#f6f6f6",
    x = 0.5, y = 0.5, width = 2, height = 1,
    just = c("left", "top"), default.units = "inches"
)

## Annotate genome label
bb_annoGenomeLabel(
    plot = genesPlot, scale = "Mb",
    x = 0.5, y = 1.5,
    just = c("left", "top"), default.units = "inches"
)

## ----philosphy_6, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE----------
## Load data
data("bb_imrHicData")
data("bb_bedpeData")
data("bb_imrH3K27acData")

## Define genomic region and widths of all plots
params <- bb_params(
    chrom = "chr21", chromstart = 28000000, chromend = 30300000,
    x = unit(0.5, "inches"), width = unit(2, "inches")
)

## Create BentoBox page
bb_pageCreate(width = 3, height = 4, default.units = "inches")

## Plot Hi-C data in region
hicPlot <- bb_plotHicSquare(
    data = bb_imrHicData,
    params = params,
    y = 0.5, height = 2,
    just = c("left", "top"), default.units = "inches"
)


## Plot and align DNA loops in region
bedpeLoops <- bb_plotPairsArches(
    data = bb_bedpeData,
    params = params, fill = "black", linecolor = "black",
    y = 2.6, height = 0.45,
    just = c("left", "top"), default.units = "inches"
)

## Plot and align signal track in region
signalPlot <- bb_plotSignal(
    data = bb_imrH3K27acData,
    params = params,
    y = 3, height = 0.45,
    just = c("left", "top"), default.units = "inches"
)

## Annotate genome label
bb_annoGenomeLabel(
    plot = signalPlot,
    params = params, scale = "Mb",
    y = 3.47, just = c("left", "top"), default.units = "inches"
)

## ----eval=TRUE, echo = TRUE, message=FALSE,warning=FALSE----------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ggplot2)
data(mtcars)
data("bb_imrHicData")
head(bb_imrHicData)

## ----eval=FALSE---------------------------------------------------------------
#  ## Create a BentoBox page
#  bb_pageCreate(width = 7, height = 3.5, default.units = "inches")
#  
#  ## Define genomic region with `bb_params`
#  params <- bb_params(chrom = "chr21",
#                      chromstart = 28000000, chromend = 30300000)
#  
#  ## Plot Hi-C data
#  hicPlot <- bb_plotHicTriangle(
#      data = bb_imrHicData,
#      params = params,
#      x = 2, y = 0.5, width = 3, height = 2,
#      just = "top", default.units = "inches"
#  )
#  
#  ## Add color scale annotation
#  bb_annoHeatmapLegend(
#      plot = hicPlot,
#      x = 3.5, y = 0.5, width = 0.12, height = 1,
#      just = c("right", "top"), default.units = "inches"
#  )
#  
#  ## Plot gene track in same genomic region
#  genesPlot <- bb_plotGenes(
#      params = params,
#      x = 0.5, y = 2.25, width = 3, height = 0.5,
#      just = c("left", "top"), default.units = "inches"
#  )
#  
#  ## Label genomic region
#  bb_annoGenomeLabel(
#      plot = genesPlot,
#      x = 0.5, y = 2.8, just = c("left", "top"), default.units = "inches"
#  )
#  
#  ## Create and place mtcars boxplot
#  boxPlot <- ggplot(mtcars) +
#      geom_boxplot(aes(gear, disp, group = gear))
#  bb_plotGG(
#      plot = boxPlot,
#      x = 4, y = 0.5, width = 2.5, height = 2,
#      just = c("left", "top"), default.units = "inches"
#  )
#  
#  ## Label panels
#  bb_plotText(
#      label = "A", fontsize = 16, fontface = "bold",
#      x = 0.75, y = 0.5, just = "center", default.units = "inches"
#  )
#  bb_plotText(
#      label = "B", fontsize = 16, fontface = "bold",
#      x = 3.75, y = 0.5, just = "center", default.units = "inches"
#  )
#  
#  ## Hide page guides
#  bb_pageGuideHide()

## ----philosophy_7, eval=TRUE, message=FALSE, echo=FALSE, fig.width=8----------
## Create a BentoBox page
bb_pageCreate(
    width = 7, height = 3.5, default.units = "inches",
    xgrid = 0, ygrid = 0, showGuides = FALSE
)

## Define genomic region with `bb_params`
params <- bb_params(chrom = "chr21", 
                    chromstart = 28000000, chromend = 30300000)

## Plot Hi-C data
hicPlot <- bb_plotHicTriangle(
    data = bb_imrHicData,
    params = params,
    x = 2, y = 0.5, width = 3, height = 2,
    just = "top", default.units = "inches"
)

## Add color scale annotation
bb_annoHeatmapLegend(
    plot = hicPlot,
    x = 3.5, y = 0.5, width = 0.12, height = 1,
    just = c("right", "top"), default.units = "inches"
)

## Plot gene track in same genomic region
genesPlot <- bb_plotGenes(
    params = params,
    x = 0.5, y = 2.25, width = 3, height = 0.5,
    just = c("left", "top"), default.units = "inches"
)

## Label genomic region
bb_annoGenomeLabel(
    plot = genesPlot,
    x = 0.5, y = 2.8, just = c("left", "top"), default.units = "inches"
)

## Create and place mtcars boxplot
boxPlot <- ggplot(mtcars) +
    geom_boxplot(aes(gear, disp, group = gear))
bb_plotGG(
    plot = boxPlot,
    x = 4, y = 0.5, width = 2.5, height = 2,
    just = c("left", "top"), default.units = "inches"
)

## Label panels
bb_plotText(
    label = "A", fontsize = 16, fontface = "bold",
    x = 0.75, y = 0.5, just = "center", default.units = "inches"
)
bb_plotText(
    label = "B", fontsize = 16, fontface = "bold",
    x = 3.75, y = 0.5, just = "center", default.units = "inches"
)

## ----eval=FALSE---------------------------------------------------------------
#  ## Change colors while maintaining plot layout
#  hicPlot <- bb_plotHicTriangle(
#      data = bb_imrHicData,
#      params = params,
#      palette = colorRampPalette(c("white", "steel blue")),
#      x = 2, y = 0.5, width = 3, height = 2,
#      just = "top", default.units = "inches"
#  )
#  
#  genesPlot <- bb_plotGenes(
#      params = params,
#      fill = c("grey", "grey"), fontcolor = c("black", "black"),
#      x = 0.5, y = 2.25, width = 3, height = 0.5,
#      just = c("left", "top"), default.units = "inches"
#  )

## ----philosophy_8, eval=TRUE, message=FALSE, echo=FALSE, fig.width=8----------
## Create a BentoBox page
bb_pageCreate(
    width = 7, height = 3.5, default.units = "inches",
    xgrid = 0, ygrid = 0, showGuides = FALSE
)

## Define genomic region with `bb_params`
params <- bb_params(chrom = "chr21", 
                    chromstart = 28000000, chromend = 30300000)

## Plot Hi-C data
hicPlot <- bb_plotHicTriangle(
    data = bb_imrHicData,
    params = params,
    palette = colorRampPalette(c("white", "steel blue")),
    x = 2, y = 0.5, width = 3, height = 2,
    just = "top", default.units = "inches"
)

## Add color scale annotation
bb_annoHeatmapLegend(
    plot = hicPlot,
    x = 3.5, y = 0.5, width = 0.12, height = 1,
    just = c("right", "top"), default.units = "inches"
)

## Plot gene track in same genomic region
genesPlot <- bb_plotGenes(
    params = params,
    fill = c("grey", "grey"), fontcolor = c("black", "black"),
    x = 0.5, y = 2.25, width = 3, height = 0.5,
    just = c("left", "top"), default.units = "inches"
)

## Label genomic region
bb_annoGenomeLabel(
    plot = genesPlot,
    x = 0.5, y = 2.8, just = c("left", "top"), default.units = "inches"
)

## Create and place mtcars boxplot
boxPlot <- ggplot(mtcars) +
    geom_boxplot(aes(gear, disp, group = gear))
bb_plotGG(
    plot = boxPlot,
    x = 4, y = 0.5, width = 2.5, height = 2,
    just = c("left", "top"), default.units = "inches"
)

## Label panels
bb_plotText(
    label = "A", fontsize = 16, fontface = "bold",
    x = 0.75, y = 0.5, just = "center", default.units = "inches"
)
bb_plotText(
    label = "B", fontsize = 16, fontface = "bold",
    x = 3.75, y = 0.5, just = "center", default.units = "inches"
)

