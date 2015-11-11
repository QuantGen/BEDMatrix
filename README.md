BEDMatrix
=========

[![Travis-CI Build Status](https://travis-ci.org/QuantGen/BEDMatrix.svg?branch=master)](https://travis-ci.org/QuantGen/BEDMatrix)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BEDMatrix)](http://cran.r-project.org/package=BEDMatrix)

BEDMatrix is an R package that provides a wrapper around [binary PED (also known as BED) files](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed), one of the genotype/phenotype file formats of [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/), the whole genome association analysis toolset. BEDMatrix objects are created in R by simply providing the path to a BED file and once created, they behave similarly to regular matrices with the advantage that genotypes are retrieved on demand without loading the entire file into memory. This allows handling of very large files with limited use of memory. Technically, a BEDMatrix is a memory-mapped matrix backed by a binary PED file.


Example
-------

This example uses a very simple BED file that is bundled with the R package. It was generated from the PLINK files in the [`inst/extdata` folder](https://github.com/QuantGen/BEDMatrix/tree/master/inst/extdata) ([see below on how to convert a PED file to a BED file](#how-do-i-create-a-bed-file-from-a-ped-file-using-plink)).

```r
# Get path to example BED file
path <- system.file('extdata', 'example.bed', package = 'BEDMatrix')

# Wrap example BED file in matrix
m <- BEDMatrix(path)

# Print matrix
m[]

# Print dimensions
dim(m)
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


Usage
-----

1. Load the library: `library(BEDMatrix)`
2. Create a new BEDMatrix object by passing in the path to the binary PED file (`path`): `m <- BEDMatrix('plink.bed')`. The `BEDMatrix` constructor will try to parse a FAM and BIM file of the same name as the BED file to determine the number and names of individuals and the number and names of SNPs. If either one of those files are not present, it is necessary to provide the number of individuals (`n`) and the number of SNPs (`p`) explicitly as parameters of the function: `m <- BEDMatrix(path = 'plink.bed', n = 100, p = 10000)` Passing `n` and `p` manually is also helpful if the FAM and BIM files are large and parsing them would take too long.
3. Extract information from the BEDMatrix as if it were a regular matrix, e.g. `m[, 3]`
4. Report any missing functionality or bugs: https://github.com/QuantGen/BEDMatrix/issues/new :)


FAQ
---

### How do I create a BED file from a PED file using PLINK?

```
plink --file myfile --make-bed
```

### Creating BEDMatrix objects is slow. How can I speed up the process?

Pass `n` and `p` manually. That way, the BIM and FAM files are not parsed.
