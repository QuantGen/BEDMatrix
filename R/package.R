#' BEDMatrix: A Wrapper for Binary PED Files
#' 
#' The BEDMatrix package provides a wrapper around binary PED files (so-called
#' BED files) that behaves just like a regular R matrix, but retrieves
#' genotypes on demand without loading the entire BED file into memory. The
#' goal is to support huge PED files, and to save time on initially reading in
#' PED files into the R environment.
#' 
#' @docType package
#' @name BEDMatrix-package
#' @useDynLib BEDMatrix
#' @import methods Rcpp
#' @aliases NULL
NULL


.onAttach <- function(libname, pkgname) {
    packageStartupMessage("The BEDMatrix package was supported by the National Institutes of Health (Grant: R01GM101219, R01GM099992).")
    packageStartupMessage()
}


loadModule("mod_BEDMatrix", TRUE)


release_questions <- function() {
    c("Have you updated the NEWS file?")
}
