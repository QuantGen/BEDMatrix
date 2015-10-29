BEDMatrix
=========

[![Travis-CI Build Status](https://travis-ci.org/QuantGen/BEDMatrix.svg?branch=master)](https://travis-ci.org/QuantGen/BEDMatrix)

The BEDMatrix package provides a wrapper around [binary PED files](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed) (so-called BED files) that behaves just like a regular R matrix, but retrieves genotypes on demand without loading the entire BED file into memory. The goal is to support huge PED files, and to save time on initially reading in PED files into the R environment.


Example
-------

This example uses a very simple BED file that is bundled with the R package. It was generated from the PLINK files in the [`inst/extdata` folder](https://github.com/QuantGen/BEDMatrix/tree/master/inst/extdata).

```r
# Get path to example BED file
path <- system.file('extdata', 'example.bed', package = 'BEDMatrix')

# Wrap example BED file in matrix
m <- BEDMatrix(path)

# Print matrix
m[]
```


Installation
------------

BEDMatrix depends on [Boost.IOStreams](http://www.boost.org/doc/libs/1_59_0/libs/iostreams/doc/index.html) for memory-mapped files. This is not a header-only library and is therefore not part of the [BH](https://cran.r-project.org/web/packages/BH/) package for R. You will need to install Boost on your system, for example by running `sudo apt install libboost-all-dev` on Ubuntu or `brew install boost` on OS X. Windows is currently unsupported.

The package is not available on CRAN yet. In the meantime, it can be installed using the [devtools](https://github.com/hadley/devtools) package directly from GitHub.

1. Install devtools: `install.packages('devtools')`
2. Load devtools: `library(devtools)`
3. Download BEDMatrix: `install_github('QuantGen/BEDMatrix')`


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

### There are errors when installing the package. What do I do now?

Did you install Boost? See [Installation](#installation) section above. If you are on Windows, see [the question below](#how-do-i-install-this-package-on-windows). If Boost is installed and it still doesn't work, please [create an issue](https://github.com/QuantGen/BEDMatrix/issues/new) and let me know.

### How do I install this package on Windows?

I don't know, but please get in touch if you find out!
