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
