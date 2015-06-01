# BEDMatrix
A wrapper for binary PED files.

## Installation
1. Install devtools: `install.packages('devtools')`
2. Load devtools: `library(devtools)`
3. Download BEDMatrix: `install_github('QuantGen/BEDMatrix')`

## Usage
1. Load the library: `library(BEDMatrix)`
2. Create a new BEDMatrix object by passing in a path to the binary PED file (`path`), the number of individuals (`n`), and the number of SNPs (`p`): `m <- BEDMatrix(path = 'plink.bed', n = 100, p = 10000)`
3. Extract information from the BEDMatrix as if it were a regular matrix, e.g. `m[, 3]`
4. Report any missing functionality or bugs: https://github.com/QuantGen/BEDMatrix/issues/new :)
