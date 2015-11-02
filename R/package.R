#' @useDynLib BEDMatrix
#' @import methods Rcpp
NULL


.onAttach <- function(libname, pkgname) {
    packageStartupMessage("The BEDMatrix package was supported by the National Institutes of Health (Grant: R01GM101219, R01GM099992).")
}


loadModule("mod_BEDMatrix", TRUE)
