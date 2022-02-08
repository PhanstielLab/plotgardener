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
