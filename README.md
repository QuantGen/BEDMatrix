BEDMatrix
=========

[![Travis-CI Build Status](https://travis-ci.org/QuantGen/BEDMatrix.svg?branch=master)](https://travis-ci.org/QuantGen/BEDMatrix)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BEDMatrix)](http://cran.r-project.org/package=BEDMatrix)

BEDMatrix is an R package that provides a wrapper around [binary PED (also known as BED) files](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed), one of the genotype/phenotype file formats of [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/), the whole genome association analysis toolset. BEDMatrix objects are created in R by simply providing the path to a BED file and once created, they behave similarly to regular matrices with the advantage that genotypes are retrieved on demand without loading the entire file into memory. This allows handling of very large files with limited use of memory. Technically, a BEDMatrix is a memory-mapped matrix backed by a binary PED file.

This package is deliberately kept simple. For computational methods that use BEDMatrix check out the [BGData package](https://github.com/QuantGen/BGData).


Example
-------

This example uses a BED version of the `mice` dataset that is included in the [BGLR package](https://github.com/gdlc/BGLR-R). See the [mice.bed gist](https://gist.github.com/agrueneberg/812564cbe860db4ee864d019be940aaf) for instructions on how it was generated.

To get the path to the example BED file (`system.file` finds the full file names of files in packages and is only used to find the example data):

```r
> path <- system.file("extdata", "mice.bed", package = "BEDMatrix")
```

To wrap the example BED file in a BEDMatrix object:

```r
> m <- BEDMatrix(path)
Extracting number of individuals and rownames from FAM file...
Extracting number of markers and colnames from BIM file...
```

To get the dimensions of the BEDMatrix object:

```r
> dim(m)
[1]  1814 10346
```

To extract a subset of the BEDMatrix object:

```r
> m[1:10, 1:3]
                      rs3683945_G rs3707673_G rs6269442_G
A048005080_A048005080           1           1           1
A048006063_A048006063           1           1           2
A048006555_A048006555           2           0           2
A048007096_A048007096           1           1           1
A048010273_A048010273           2           0           2
A048010371_A048010371           1           1           1
A048011040_A048011040           1           1           1
A048011287_A048011287           1           1           2
A048011567_A048011567           1           1           2
A048013559_A048013559           2           0           2
```


Installation
------------

To get the current released version from CRAN:

```r
install.packages("BEDMatrix")
```

To get the current development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("QuantGen/BEDMatrix")
```


Documentation
-------------

Further documentation can be found on [RDocumentation](http://www.rdocumentation.org/packages/BEDMatrix).
