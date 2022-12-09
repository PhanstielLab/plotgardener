# CHANGES IN VERSION 1.4.2
## BUG FIXES

* `annoPixels` detects and annotates all pixels for `plotHicRectangle` plots.

# CHANGES IN VERSION 1.4.1
## BUG FIXES

* yscales for `plotHicRectangle` and `plotHicTriangle` reflect distance in 
Hi-C bins.

# VERSION 1.4.0

Version bump for Bioconductor 3.16 release.

# CHANGES IN VERSION 1.3.8 AND 1.3.9
## BUG FIXES

* Fixed vignette links in "Introduction to plotgardener" vignette.

# CHANGES IN VERSION 1.3.7
## BUG FIXES 

* `plotGenes` and related functions will appropriately check for and handle
custom `OrgDb`s.
* `getExons` will double-check for appropriate chromosome data to avoid
incorrect plotting based on related chromosome contigs.

## NEW FEATURES

* `plotManhattan` p-value data can be scaled according to a custom 
transformation, rather than being limited to -log10.
* `plotRanges` elements can be ordered randomly or by decreasing width before
plotted row assignment.

# CHANGES IN VERSION 1.3.6
## BUG FIXES 

* `mapColors` can appropriately map colors to a numeric vector with the same 
values, so long as a range is provided.

# CHANGES IN VERSION 1.3.5
## NEW FEATURES

* `plotMultiSignal` function can plot multiple signal track data sets in line
with each other.
* `calcSignalRange` helper function will calculate an appropriate range for
multiple signal data sets.
* `pageLayoutRow` and `pageLayoutCol` functions for generating row and column
positions for a number of plot elements.
* Gene transcripts can be highlighted by gene name or transcript name with the
parameter `transcriptHighlights` in `plotTranscripts`.

# CHANGES IN VERSION 1.3.4
## BUG FIXES

* `plotPairsArches` Bezier curve height calculations were fixed for pairs
with different sized anchors.

# CHANGES IN VERSION 1.3.3
## NEW FEATURES

* `plotSignal` can now plot negative signal data alone or listed as a second 
file.
* A `label` parameter has been added for `plotSignal` for convenient labeling. 

# CHANGES IN VERSION 1.3.2
## BUG FIXES 

* `plotSignal` range parsing bug fixes were resolved.

* Note about double page rendering has been added to `pageGuideHide()`
documentation.

# CHANGES IN VERSION 1.3.1

* `plotSignal` bug fixes related to function not finding posSignal2 and 
negSignal2 variables with insufficient data.

* Documentation to introduction vignette has been added to explain double 
page rendering when using any removal function, particularly `pageGuideHide()`.


# VERSION 1.3.0

Version bump for Bioconductor 3.15 release.

# CHANGES IN VERSION 1.1.18
## NEW FEATURES

* `hicTriangles` and `hicRectangles` can now be annotated with `annoDomains`
or `annoPixels` if they are flipped.


# CHANGES IN VERSION 1.1.17
## NEW FEATURES

* `plotIdeogram` can now accept custom colors with a `fill` parameter. Colors
can be specified with a named or unnamed vector. To see which stains are being
assigned which colors, look inside the `ideogram` object.

# CHANGES IN VERSION 1.1.16

## BUG FIXES 

* ENTREZ IDs obtained from `AnnotationDbi::select()` are subset just for
`ENTREZID` column when determining default gene priorities, eliminating `dplyr`
incompatible types error.
* All plus and minus strand gene name label parsing in `plotGenes` is now 
carried out only if there is a non-zero number of that strand's 
genes.

# CHANGES IN VERSION 1.1.15

* Citation linked for `plotgardener` publication in Bioinformatics.

# CHANGES IN VERSION 1.1.14

## BUG FIXES 

* `plotSignal` yrange parsing for negative scores now has fixed the typo on
line 418 from "score2" to "score".

# CHANGES IN VERSION 1.1.13

## BUG FIXES 

* `plotSignal` default yrange parsing now catches the invalid 0,0 range and no
longer throws a viewport related error.

# CHANGES IN VERSION 1.1.12

## BUG FIXES 

* `readHic` and functions related to the reading of .hic files now leaves the
chromosome input formatted as is (e.g. "chr1" and "1"). Functions will throw an 
error if the input chromosome is not found in the chromosomes listed in 
the .hic file.

# CHANGES IN VERSION 1.1.11

## BUG FIXES 

* `annoDomains` coordinates fixed for `plotHicRectangle`.
* Clipping logic for `plotPairsArches` now clips arches both on left and right
side of plot.
* Subsetting `plotPairs` logic fixed to match `plotPairsArches`.
    
## NEW FEATURES

* `clip.noAnchors` parameter in `plotPairsArches` allows for inclusion or 
clipping of arches that do not have anchors in the given genomic region.
* `plotPairsArches` now allows for column name input to designate `archHeights`.


# CHANGES IN VERSION 0.99.11

## BUG FIXES 

* `:` added back to readHic for `strawr` region parsing.
* `plotHicSquare` subsetting fixed for off diagonal regions.
    
## NEW FEATURES

* All Hi-C functions now allow input of remote Hi-C files.

# CHANGES IN VERSION 0.99.9

This package was previously called `BentoBox`.

# CHANGES IN VERSION 0.99.0

## NEW FEATURES

* Version bump to 0.99.0 for Bioconductor package submission.
* `bb_mapColors` function for users to map a vector to a palette of colors.
* `linecolor` parameter in `bb_plotPairs`, `bb_plotPairsArches`,
and `bb_plotRanges` now accepts a single value, a vector of colors, 
a `colorby` object, or the value "fill".

# CHANGES IN VERSION 0.0.0.14

## BUG FIXES 

* R version requirement changed to (R >= 4.1.0) for proper plot placement.
    
## NEW FEATURES

* `colorby` object now has a `scalePerRegion` parameter to scale numerical
`colorby` data to the range of data in a plotted genomic region.

# CHANGES IN VERSION 0.0.0.13

## SIGNIFICANT USER-VISIBLE CHANGES

* `bb_plotManhattan` `fill` paramete now accepts a single value, 
a vector of colors, or a `colorby` object.

# CHANGES IN VERSION 0.0.0.12

## SIGNIFICANT USER-VISIBLE CHANGES

* `colorby` constructor now includes optional palette specification.
* `bb_plotPairs`, `bb_plotPairsArches`, and `bb_plotRanges` `fill` parameter
now accepts a single value, a vector of colors, or a `colorby` object.
    
## BUG FIXES    

* `GInteractions` assembly match checking moved before dataframe 
conversion.

# CHANGES IN VERSION 0.0.0.11

## SIGNIFICANT USER-VISIBLE CHANGES

* Data moved to `plotgardenerData` package.
* Default genome assembly updated to "hg38".
    
## BUG FIXES

* Streamlined parameter parsing and data reading logic.

# CHANGES IN VERSION 0.0.0.10

## NEW FEATURES

* Added unit tests with `testthat`.
* `bb_annoDomains` function addition.
* `bb_plotSignal` vertical orientation.
    
# CHANGES IN VERSION 0.0.0.9

## NEW FEATURES

* Added a `NEWS` file to track changes to the package.
    
## BUG FIXES

* Updated viewport parsing for package `grid` version 4.1.0.
