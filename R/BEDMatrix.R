#' @useDynLib BEDMatrix
#' @importFrom Rcpp sourceCpp
NULL

#' @export
BEDMatrix <- function (path, n, p) {
  if (!file.exists(path)) {
    stop('File not found.');
  }
  n <- as.integer(n)
  p <- as.integer(p)
  obj <- list()
  class(obj) <- 'BEDMatrix'
  attr(obj, 'path') <- path
  attr(obj, 'n') <- n
  attr(obj, 'p') <- p
  attr(obj, 'dnames') <- list(
    paste0('id_', 1:n),
    paste0('mrk_', 1:p)
  )
  return(obj)
}

#' @export
print.BEDMatrix <- function (x, ...) {
  n <- attr(x, 'n')
  p <- attr(x, 'p')
  cat(paste(n, 'x', p, 'BEDMatrix'))
}

#' @export
`[.BEDMatrix` <- function (x, i, j, drop) {
  path <- attr(x, 'path')
  n <- attr(x, 'n')
  p <- attr(x, 'p')
  if (missing(i)) {
    i <- 1:n
  }
  if (missing(j)) {
    j <- 1:p
  }
  subset <- subsetBED(x, i, j)
  return(subset)
}

#' @export
dim.BEDMatrix <- function (x) {
  n <- attr(x, 'n')
  p <- attr(x, 'p')
  return(c(n, p))
}

#' @export
dimnames.BEDMatrix <- function (x) {
  dnames <- attr(x, 'dnames')
  return(dnames)
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
  attr(x, 'dnames') <- lapply(value, function (v) {
    if (!is.null(v)) {
      as.character(v)
    }
  })
  return(x)
}
