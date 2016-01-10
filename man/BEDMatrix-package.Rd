% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{package}
\name{BEDMatrix-package}
\title{\code{BEDMatrix}: A Wrapper for Binary PED Files}
\description{
The BEDMatrix package provides a wrapper around binary PED (also known as
BED) files, one of the genotype/phenotype file formats of
\href{http://pngu.mgh.harvard.edu/~purcell/plink/}{PLINK}, the whole genome
association analysis toolset. \code{BEDMatrix} objects are created by simply
providing the path to a BED file and once created, they behave similarly to
regular matrices with the advantage that genotypes are retrieved on demand
without loading the entire file into memory. This allows handling of very
large files with limited use of memory. Technically, a \code{BEDMatrix} is a
memory-mapped matrix backed by a binary PED file.
}
