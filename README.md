BEDMatrix
=========

[![Travis-CI Build Status](https://travis-ci.org/QuantGen/BEDMatrix.svg?branch=master)](https://travis-ci.org/QuantGen/BEDMatrix)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BEDMatrix)](https://cran.r-project.org/package=BEDMatrix)
[![Rdoc](http://www.rdocumentation.org/badges/version/BEDMatrix)](http://www.rdocumentation.org/packages/BEDMatrix)

BEDMatrix is an R package that provides a matrix-like wrapper around [.bed](https://www.cog-genomics.org/plink2/formats#bed), one of the genotype/phenotype file formats of [PLINK](https://www.cog-genomics.org/plink2), the whole genome association analysis toolset. BEDMatrix objects are created in R by simply providing the path to a .bed file and once created, they behave similarly to regular matrices with the advantage that genotypes are retrieved on demand without loading the entire file into memory. This allows handling of very large files with limited use of memory.

This package is deliberately kept simple. For computational methods that use BEDMatrix check out the [BGData package](https://github.com/QuantGen/BGData).


Example
-------

This example uses a dummy .bed file that is bundled with this R package. It was generated using `plink --dummy 500 1000 0.02 acgt --seed 4711 --out example` with PLINK 1.90 beta 3.452.

To get the path to the example .bed file (`system.file` finds the full file names of files in packages and is only used to find the example data):

```R
> path <- system.file("extdata", "example.bed", package = "BEDMatrix")
```

To wrap the example .bed file in a BEDMatrix object:

```R
> m <- BEDMatrix(path)
Extracting number of samples and rownames from FAM file...
Extracting number of variants and colnames from BIM file...
```

To get the dimensions of the BEDMatrix object:

```R
> dim(m)
[1] 50 1000
```

To extract a subset of the BEDMatrix object:

```R
> m[1:3, 1:5]
          snp0_A snp1_C snp2_G snp3_G snp4_G
per0_per0      0      1      1      1      0
per1_per1      1      1      1      1     NA
per2_per2      1      0      0      2      0
```


Installation
------------

Install the stable version from CRAN:

```R
install.packages("BEDMatrix")
```

Alternatively, install the development version from GitHub:

```R
# install.packages("devtools")
devtools::install_github("QuantGen/BEDMatrix")
```


Contribute
----------

- Issue Tracker: https://github.com/QuantGen/BEDMatrix/issues
- Source Code: https://github.com/QuantGen/BEDMatrix


Documentation
-------------

Further documentation can be found on [RDocumentation](http://www.rdocumentation.org/packages/BEDMatrix).
