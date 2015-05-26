#' @useDynLib BEDMatrix
#' @importFrom Rcpp sourceCpp
NULL

#' @export
setClass('BEDMatrix', slots = list(path = 'character', n = 'integer', p = 'integer', dimnames = 'list'))

#' @export
setMethod('initialize', signature(.Object = 'BEDMatrix'), function (.Object, path, n, p) {
  .Object@path <- path
  .Object@n <- as.integer(n)
  .Object@p <- as.integer(p)
  .Object@dimnames <- list(NULL, NULL)
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

#' @export
dim.BEDMatrix <- function (x) {
  return(c(x@n, x@p))
}

#' @export
dimnames.BEDMatrix <- function (x) {
  return(x@dimnames)
}

#' @export
`dimnames<-.BEDMatrix` <- function (x, value) {
  d <- dim(x)
  v1 <- value[[1]]
  v2 <- value[[2]]
  if (!is.list(value) || length(value) != 2 ||
      !(is.null(v1) || length(v1) == d[1]) ||
      !(is.null(v2) || length(v2) == d[2])) {
    stop('invalid dimnames')
  }
  x@dimnames <- lapply(value, function (v) {
    if (!is.null(v)) {
      as.character(v)
    }
  })
  return(x)
}
