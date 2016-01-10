BEDMatrix
=========

[![Travis-CI Build Status](https://travis-ci.org/QuantGen/BEDMatrix.svg?branch=master)](https://travis-ci.org/QuantGen/BEDMatrix)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BEDMatrix)](http://cran.r-project.org/package=BEDMatrix)

BEDMatrix is an R package that provides a wrapper around [binary PED (also known as BED) files](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed), one of the genotype/phenotype file formats of [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/), the whole genome association analysis toolset. BEDMatrix objects are created in R by simply providing the path to a BED file and once created, they behave similarly to regular matrices with the advantage that genotypes are retrieved on demand without loading the entire file into memory. This allows handling of very large files with limited use of memory. Technically, a BEDMatrix is a memory-mapped matrix backed by a binary PED file.


Example
-------

This example uses a very simple BED file that is bundled with the R package. It was generated from the PLINK files in the [`inst/extdata` folder](https://github.com/QuantGen/BEDMatrix/tree/master/inst/extdata) ([see below on how to convert a PED file to a BED file](#how-do-i-create-a-bed-file-from-a-ped-file-using-plink)).

To get the path to the example BED file (`system.file` finds the full file names of files in packages and is only used to find the example data):

```r
> path <- system.file("extdata", "example.bed", package = "BEDMatrix")
```

To wrap the example BED file in a BEDMatrix object:

```r
> m <- BEDMatrix(path)
Extracting number of individuals and rownames from FAM file...
Extracting number of markers and colnames from BIM file...
```

To print the dimensions of the BEDMatrix object:

```r
> dim(m)
[1] 6 3
```

To print the entire BEDMatrix object:

```r
> m[]
    snp1_G snp2_1 snp3_A
1_1      2      0      0
1_2      0     NA      1
1_3     NA      1      1
2_1      0      0     NA
2_2      0      0     NA
2_3      0      0      2
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
