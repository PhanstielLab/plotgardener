
# <img src="man/figures/bento-logotype-singleline-black.png" width="400px" style="background-color:white;border=transparent" /> <img src="man/figures/logo.png" align="right" width="140px" style="padding-left:20px; background-color:white" />

<!-- badges: start -->
[![R build status](https://github.com/PhanstielLab/BentoBox/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/PhanstielLab/BentoBox/actions)
<!-- badges: end -->

## Overview

`BentoBox` is a genomic data visualization package for R. Using `grid`
graphics, `BentoBox` empowers users to programmatically and flexibly
generate multi-panel figures. `BentoBox` accomplishes these goals by
utilizing 1) a coordinate-based plotting system, and 2) edge-to-edge
containerized data visualization. The coordinate-based plotting system
grants users precise control over the size, position, and arrangement of
plots. Its edge-to-edge plotting functions preserve the mapping between
user-specified containers and the represented data. This allows users to
stack plots with confidence that vertically aligned data will correspond
to the same regions. For more information about BentoBox’s philosophy
and design, check out the `Our Philosophy` page.

Specialized for genomic data, `BentoBox` also contains functions to read
and plot multi-omic data quickly and easily. These functions are
integrated with Bioconductor packages to flexibly accommodate a large variety 
of genomic assemblies. `BentoBox` can address an
endless number of use cases, including: dynamic exploration of genomic
data, arrangment into multi-omic layouts, and survey plotting for
quickly viewing data across the genome. Check out our `vignettes` for
detailed examples and suggested use cases!

## Installation

BentoBox can be installed from GitHub as follows:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("remotes", quietly = TRUE))
    BiocManager::install("remotes")

remotes::install_github("PhanstielLab/BentoBox")
packageVersion("BentoBox")
```

Cytoband annotation datasets and example datasets and files are included
with the package BentoBoxData:

``` r
remotes::install_github("PhanstielLab/BentoBoxData")
```

## Usage

<img src="man/figures/homePage-1.png" width="672" style="display: block; margin: auto;" />

``` r
## Load libraries and datasets
library("BentoBox")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("BentoBoxData")
library("AnnotationHub")
data("GM12878_HiC_10kb")
data("IMR90_HiC_10kb")
data("GM12878_ChIP_CTCF_signal")
data("IMR90_ChIP_CTCF_signal")
data("GM12878_ChIP_H3K27ac_signal")
data("IMR90_ChIP_H3K27ac_signal")

## Create a BentoBox page
bbPageCreate(width = 7, height = 4.25, default.units = "inches")

##########################
######## Panel A #########
##########################
## Text section label
bbPlotText(label = "A", fontsize = 12,
            x = 0.25, y = 0.25, just = "left", default.units = "inches")
## Set genomic and dimension parameters in a `bbParams` object
params_a <- bbParams(chrom = "chr21", chromstart = 28000000, chromend = 30300000, 
                      assembly = "hg19",
                      x = 0.25, width = 2.75, just = c("left", "top"), default.units = "inches")
## Double-sided Hi-C Plot
hicPlot_top <- bbPlotHicSquare(data = GM12878_HiC_10kb, params = params_a,
                                zrange = c(0, 200), resolution = 10000,
                                half = "top",
                                y = 0.5, height = 2.75)
hicPlot_bottom <- bbPlotHicSquare(data = IMR90_HiC_10kb, params = params_a,
                                   zrange = c(0, 70), resolution = 10000,
                                   half = "bottom",
                                   y = 0.5, height = 2.75)
## Annotate Hi-C heatmap legends
bbAnnoHeatmapLegend(plot = hicPlot_bottom, fontsize = 7, 
                     x = 3.05, y = 0.5,
                     width = 0.07, height = 0.5, just = c("left", "top"),
                     default.units = "inches")
bbAnnoHeatmapLegend(plot = hicPlot_top, fontsize = 7, 
                     x = .125, y = 0.5,
                     width = 0.07, height = 0.5, just = c("left", "top"),
                     default.units = "inches")
## Plot gene track
genes_a <- bbPlotGenes(params = params_a, stroke = 1, fontsize = 6,
                        y = 3.35, height = 0.4)
## Annotate genome label
bbAnnoGenomeLabel(plot = genes_a, params = params_a, 
                   scale = "Mb", fontsize = 7,
                   y = 3.85)
##########################
######## Panel B #########
##########################
## Text section label
bbPlotText(label = "B", fontsize = 12,
            x = 3.5, y = 0.25, just = "left", default.units = "inches")
## Plot ideogram
bbPlotIdeogram(chrom = "chr21", assembly = "hg19",
                x = 3.5, y = 0.5,
                width = 3.25, height = 0.15, just = c("left", "top"), default.units = "inches")
## Add text to ideogram
bbPlotText(label = "Chromosome 21", fontsize = 8, fontcolor = "darkgrey",
            x = 6.75, y = 0.4, just = "right", default.units = "inches")
##########################
######## Panel C #########
##########################
## Text section label
bbPlotText(label = "C", fontsize = 12,
            x = 3.5, y = 1, just = c("left", "top"), default.units = "inches")
## Set genomic and dimension parameters in a `bbParams` object
params_c <- bbParams(chrom = "chr21", chromstart = 28150000, chromend = 29150000, 
                      assembly = "hg19",
                      x = 3.5, width = 1.5, default.units = "inches")
## Set signal track data ranges
ctcf_range <- bbParams(range = c(0, 77),
                        assembly = "hg19")
hk_range <- bbParams(range = c(0, 32.6),
                      assembly = "hg19")
## Plot Hi-C triangle
hic_gm <- bbPlotHicTriangle(data = GM12878_HiC_10kb, params = params_c,
                             zrange = c(0, 200), resolution = 10000,
                             y = 1.75, height = 0.75, just = c("left", "bottom"))
## Annotate Hi-C heatmap legend
bbAnnoHeatmapLegend(plot = hic_gm, fontsize = 7, 
                     x = 5, y = 1, width = 0.07, height = 0.5, 
                     just = c("right", "top"), default.units = "inches")
## Plot CTCF signal
ctcf_gm <- bbPlotSignal(data = GM12878_ChIP_CTCF_signal, params = c(params_c, ctcf_range),
                         fill = "#253494", linecolor = "#253494",
                         y = 1.95, height = 0.6)
## CTCF label
bbPlotText(label = "CTCF", fontcolor = "#253494", fontsize = 8,
            x = 3.5, y = 1.95, just = c("left","top"), default.units = "inches")
## Plot H3K27ac signal
hk_gm <- bbPlotSignal(data = GGM12878_ChIP_H3K27ac_signal, params = c(params_c, hk_range),
                       fill = "#37a7db", linecolor = "#37a7db",
                       y = 3.25, height = 0.6, just = c("left", "bottom"))
## H3K27ac label
bbPlotText(label = "H3K27ac", fontcolor = "#37a7db", fontsize = 8,
            x = 3.5, y = 2.65, just = c("left","top"), default.units = "inches")
## Plot genes
genes_gm <- bbPlotGenes(params = params_c, stroke = 1, fontsize = 6,
                         strandLabels = FALSE,
                         y = 3.35, height = 0.4)
## Annotate genome label
bbAnnoGenomeLabel(plot = genes_gm, params = params_c, 
                   scale = "Kb", fontsize = 7,
                   y = 3.85)
##########################
######## Panel D #########
##########################
## Text section label
bbPlotText(label = "D", fontsize = 12,
            x = 5.25, y = 1, just = c("left", "top"), default.units = "inches")
## Set genomic and dimension parameters in a `bbParams` object
params_d <- bbParams(chrom = "chr21", chromstart = 28150000, chromend = 29150000, 
                      assembly = "hg19",
                      x = 6.75, width = 1.5, default.units = "inches")
## Plot Hi-C triangle
hic_imr <- bbPlotHicTriangle(data = IMR90_HiC_10kb, params = params_d,
                              zrange = c(0, 70), resolution = 10000,
                              y = 1.75, height = 0.75, just = c("right", "bottom"))
## Annotate Hi-C heatmap legend
bbAnnoHeatmapLegend(plot = hic_imr, fontsize = 7,
                     x = 6.75, y = 1, width = 0.07, height = 0.5, just = c("right", "top"))
## Plot CTCF signal
ctcf_imr <- bbPlotSignal(data = IMR90_ChIP_CTCF_signal, params = c(params_d, ctcf_range),
                          fill = "#253494", linecolor = "#253494",
                          y = 1.95, height = 0.6, just = c("right", "top"))
## Plot H3K27ac signal
hk_imr <- bbPlotSignal(data = IMR90_ChIP_H3K27ac_signal, params = c(params_d, hk_range),
                        fill = "#37a7db", linecolor = "#37a7db",
                        y = 3.25, height = 0.6, just = c("right", "bottom"))
## Plot gene track
genes_imr <- bbPlotGenes(params = params_d, stroke = 1, fontsize = 6,
                          strandLabels = FALSE,
                          y = 3.35, height = 0.4, just = c("right", "top"))
## Annotate genome label
bbAnnoGenomeLabel(plot = genes_imr, params = params_d, 
                   scale = "Kb", fontsize = 7,
                   y = 3.85, just = c("right", "top"))
## Hide page guides
bbPageGuideHide()
```

## A word of caution

`BentoBox` is incredibly flexible and functional. However, due to this
flexibility and like all programming packages, it may not always prevent
users from making unintentional mistakes. If plot sizes are entered
incorrectly or data is mishandled, it is possible to connect multi-omic
data incorrectly. Make sure you utilize package features that reduce
human error and increase re-usability of code to get the most milage out
of BentoBox.
