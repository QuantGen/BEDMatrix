#' BEDMatrix: A Wrapper for Binary PED Files
#'
#' The BEDMatrix package provides a wrapper around binary PED (also known as
#' BED) files, one of the genotype/phenotype file formats of
#' [PLINK](https://www.cog-genomics.org/plink2), the whole genome association
#' analysis toolset. [BEDMatrix-class] objects are created by simply providing
#' the path to a BED file and once created, they behave similarly to regular
#' matrices with the advantage that genotypes are retrieved on demand without
#' loading the entire file into memory. This allows handling of very large
#' files with limited use of memory. Technically, a [BEDMatrix-class] object is
#' a memory-mapped matrix backed by a binary PED file.
#'
#' @docType package
#' @name BEDMatrix-package
#' @aliases NULL
#' @import methods Rcpp
NULL


release_questions <- function() {
    c("Have you updated the NEWS file?")
}
