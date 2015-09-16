# BEDMatrix
A wrapper for [binary PED files](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed).

## Installation
1. Install devtools: `install.packages('devtools')`
2. Load devtools: `library(devtools)`
3. Download BEDMatrix: `install_github('QuantGen/BEDMatrix')`

## Usage
1. Load the library: `library(BEDMatrix)`
2. Create a new BEDMatrix object by passing in the path to the binary PED file (`path`): `m <- BEDMatrix('plink.bed')`. The `BEDMatrix` constructor will try to parse a FAM and MAP file of the same name as the BED file to determine the number and names of individuals and the number and names of SNPs. If either one of those files are not present, it is necessary to provide the number of individuals (`n`) and the number of SNPs (`p`) explicitly as parameters of the function: `m <- BEDMatrix(path = 'plink.bed', n = 100, p = 10000)`
3. Extract information from the BEDMatrix as if it were a regular matrix, e.g. `m[, 3]`
4. Report any missing functionality or bugs: https://github.com/QuantGen/BEDMatrix/issues/new :)

## Example
```r
m <- BEDMatrix(system.file('extdata', 'example.bed', package = 'BEDMatrix'))
m[]
```
