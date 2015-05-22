#' @useDynLib BEDMatrix
#' @importFrom Rcpp sourceCpp
NULL

#' @export
setClass('BEDMatrix', slots = list(path = 'character', n = 'integer', p = 'integer'))

#' @export
setMethod('initialize', signature(.Object = 'BEDMatrix'), function (.Object, path, n, p) {
  .Object@path <- path
  .Object@n <- as.integer(n)
  .Object@p <- as.integer(p)
  return(.Object)
})

#' @export
setMethod('[', signature(x = 'BEDMatrix'), function (x, i, j, drop) {
  if (missing(i)) {
    i <- 1:x@n
  }
  if (missing(j)) {
    j <- 1:x@p
  }
  subset <- subsetBED(x@path, x@n, x@p, i, j)
  return(subset)
})
