# DEVELOPMENT VERSION
## NEW FEATURES

* Added `pageExport()` for writing the current page to PDF, PNG, TIFF, JPEG,
  or SVG without a user `dev.off()` call.
* `pageGuideHide()` and `pagePlotRemove()` now defer removals until
  `pageExport()`/`pageRedraw()` instead of calling `grid.remove()` on the
  active device, avoiding extra PDF pages.
* `plotHicSquare()`, `plotHicRectangle()`, and `plotHicTriangle()` gain
  `raster = c("auto", "never", "always")` and `raster.max` for rasterizing
  large contact maps.
* `plotSignal()` gains `label.position`, `label.offset`, numeric baselines,
  and `baseline.lty`. Baselines are now drawn after signal/cutoff grobs so
  filled tracks cannot cover them.
* `plotMultiSignal()` now honors a user-supplied `range` and passes through
  label/baseline styling parameters.
* `annoYaxis()` and `annoXaxis()` gain `tcl`, `label.offset`,
  `axisLine.lwd/col/lty`, and `tick.lwd/col/lty` for flexible axis styling.
  Default tick length is now shorter (`tcl = 0.25`) and default tick/axis
  line width is thinner (`0.5`).
* `plotNeoGenes()` now uses IGV-style display modes. `displayMode = "auto"`
  collapses genes to solid blocks at whole-chromosome scale and expands to
  gene bodies, exon blocks, UTR/CDS heights, and strand arrows when zoomed
  in. It supports `chromstart`/`chromend` neo-coordinate windows, arbitrary
  GTF/GFF GENCODE versions, `TxDb`/`OrgDb` annotations, `geneOrder`,
  `labelClusters`, `geneBackground`, configurable colors, gene highlights,
  label collision suppression, row-overflow warnings, and `geneMode`
  layouts: `"expand"` draws all genes normally, `"focus"` fades non-priority
  genes with `focus.alpha`, and `"compact"` draws all genes on one line.
  GTF/GFF files are read with `readGFF()` by default and can use an optional
  disk cache.
* Added `getNeoSignal()` and `getNeoV4C()` for mapping bigWig, GRanges,
  data-frame, and virtual-4C signal tracks into neo-chromosome coordinates,
  with optional binning and aggregation.
* `plotNeoHic(flip = "horizontal")` reverses the x-axis for neo Hi-C plots.
* `annoNeoDomains()` and `plotNeoAssembly()` expose separate seam color,
  width, and line type, and `annoNeoFragmentLabels(label.colors=)` supports
  per-fragment label colors.
* Interchromosomal `plotHicSquare()` plots are no longer mirrored and store
  explicit x/y chromosome metadata for `annoGenomeLabel()`. Add
  `aspect = "fixed"` for equal basepairs-per-inch on unequal-length
  chromosomes.
* Added SV-impacted neo-chromosome reconstruction helpers:
  `readNeoAssembly()`, `mapNeoFeatures()`, `readNeoHic()`, `plotNeoHic()`,
  `plotNeoAssembly()`, `annoNeoDomains()`, and `annoNeoFragmentLabels()`.
  `readNeoHic()` stitches fragment-pair contacts from `.hic`, `.cool`, or
  `.mcool` files into a concatenated neo-chromosome coordinate system.
  `annoNeoDomains()` draws the solid outer neo-chromosome domain and dotted
  internal fragment seams used in the neo-triangle examples.

* `plotTranscripts()` adds `labelPriority` and centers labels on the
  in-bounds transcript portion.
* `plotGenes()` adds `labelClusters`/`labelClusterRot` for dense gene
  clusters.
* Added `plotNeoGenes()` for gene/transcript models mapped into neo
  coordinates, including strand flips, gene highlights, and selected labels. `readNeoHic()` now reports
  fragment-pair contact counts and warns when inter-fragment pairs have no
  source contacts.
* `assembly()` accepts an `EnsDb` object and converts it to a TxDb for
  plotting.
* `read_rangeData()` respects named `chrom`/`start`/`end` or `chrom`/`pos`
  columns in any order.
* `annoPixels()` handles factor chromosome columns with different level sets.
* `plotLegend()` documents the cairosvg/SVG rendering caveat.

## BUG FIXES

* `readCool()` now returns all pixels in the final bin of a region. The
previous pixel-range calculation stopped at the first pixel of the last bin.
* `readCool()` applies cooler balancing weights multiplicatively
(`counts * weight1 * weight2`) instead of dividing by them. This also applies
to `sweight`, `KR`, `VC`, and `VC_SQRT` weights.
* `readCoolNorms()` now returns the exact column names stored in the file
(for example `weight`, `sweight`, `KR`, or `obj_weight`) instead of mapping
`weight` to `BALANCE`.
* `readCool()` accepts chromosome names with or without a `chr` prefix and
matches regions that do not start exactly on a bin boundary.
* `readCool()` now reads interchromosomal mcool pixels when the requested
`chrom` comes after `altchrom` in the file, and expands both coordinates for
same-chromosome off-diagonal regions.
* `readCool()` supports `norm = "auto"`, and `plotHicSquare()`,
`plotHicRectangle()`, and `plotHicTriangle()` now use `norm = "auto"` by
default. For `.hic` files this resolves to `"KR"`; for `.(m)cool` files it
chooses an available weight.
* Hi-C plotting no longer requests negative `chromstart` values from
`.hic`/`.(m)cool` files when padding a region near the start of a chromosome.
* `.checkCool()` recognizes `.mcool` files with extra root-level groups.

# VERSION 1.11.2
## NEW FEATURES

* `plotgardener` now supports .(m)cool files for Hi-C data! Functions for 
reading and extracting features from .(m)cool files include:
    - `readCool`
    - `readCoolBpResolutions`
    - `readCoolNorms`
    - `readCoolChroms` 
  .(m)cool files can now be used with any Hi-C plotting function in place of
  .hic files or other data types.

# VERSION 1.11.1
## BUG FIXES

* `plotPairs` and `plotPairsArches` check and swap order of anchors.
* `plotTranscripts` `NULL` label checking bug fix.

# VERSION 1.11.0

Version bump for Bioconductor 3.19 release.

# VERSION 1.9.5
## NEW FEATURES

* `plotHicSquare` has added `yaxisDir` parameter for flipping the direction
of genomic coordinates on the y-axis.

# VERSION 1.9.4
## BUG FIXES

* `plotManhattan` accommodates `GRanges` data input.

# VERSION 1.9.3
## BUG FIXES

* `plotGenes` appropriately handles `tibble` input for `geneHighlights`.

# VERSION 1.9.1
## NEW FEATURES

* `plotIdeogram` has added `flip` parameter to allow for
flipping the ideogram so the end can be before/above the start.

## BUG FIXES

* `plotManhattan` pch mapping is compatible when number of data 
points is less than the number of levels of a `colorby` column.


# VERSION 1.9.0

Version bump for Bioconductor 3.18 release.

# VERSION 1.7.7
## BUG FIXES

* `params` parsing removes non-standard chromosomes when identifying gene 
transcripts


# VERSION 1.7.6
## BUG FIXES

* `colorby` logic parsing is patched. 

# VERSION 1.7.5
## NEW FEATURES

* `plotPairs`, `plotPairsArches`, and `plotRanges` allow `NA` fill.

# VERSION 1.7.4
## BUG FIXES

* Removed extra page creation with `annoYaxis` `axisLine = FALSE`.

# VERSION 1.7.3
## BUG FIXES

* Removed blank page creation with pdf() calls for all major plotting functions.

# VERSION 1.7.2
## BUG FIXES

* Fixed `plotMultiSignal` width and height parsing bug.

# CHANGES IN VERSION 1.7.1
## NEW FEATURES

* `plotManhattan` y-scales can be reversed for Miami plot-style layouts.

# VERSION 1.7.0

Version bump for Bioconductor 3.17 release.

# CHANGES IN VERSION 1.5.3
## BUG FIXES

* Fixed page viewport parsing bug fixes related to R version 4.3.0 updates.

# CHANGES IN VERSION 1.5.2
## BUG FIXES

* `annoPixels` detects and annotates all pixels for `plotHicRectangle` plots.

# CHANGES IN VERSION 1.5.1
## BUG FIXES

* yscales for `plotHicRectangle` and `plotHicTriangle` reflect distance in 
Hi-C bins.

# VERSION 1.5.0

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