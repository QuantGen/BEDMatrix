#' @useDynLib BEDMatrix
#' @importFrom Rcpp sourceCpp
NULL

#' @export
BEDMatrix <- function (path, n = NULL, p = NULL) {
  if (!file.exists(path)) {
    stop('File not found.');
  }
  dir <- substr(path, 1, nchar(path) - 4)
  # Check if FAM file exists.
  if (file.exists(paste0(dir, '.fam'))) {
    fam <- readLines(paste0(dir, '.fam'))
    # Determine n.
    n <- length(fam)
    # Determine rownames.
    rownames <- sapply(strsplit(fam, ' '), function (line) {
      return(paste0(line[1], '_', line[2]))
    })
  } else {
    if (is.null(n)) {
      stop('FAM file of same name not found. Provide number of individuals (n).')
    }
    rownames <- paste0('id_', 1:n)
  }
  # Check if MAP file exists.
  if (file.exists(paste0(dir, '.map'))) {
    map <- readLines(paste0(dir, '.map'))
    # Determine p.
    p <- length(map)
    # Determine colnames.
    colnames <- sapply(strsplit(map, ' '), function (line) {
      return(line[2])
    })
  } else {
    if (is.null(p)) {
      stop('MAP file of same name not found. Provide number of markers (p).')
    }
    colnames <- paste0('mrk_', 1:p)
  }
  n <- as.integer(n)
  p <- as.integer(p)
  obj <- list()
  class(obj) <- 'BEDMatrix'
  attr(obj, 'path') <- path
  attr(obj, 'n') <- n
  attr(obj, 'p') <- p
  attr(obj, 'dnames') <- list(
    rownames,
    colnames
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
`[.BEDMatrix` <- function (x, i, j, drop = TRUE) {
  path <- attr(x, 'path')
  n <- attr(x, 'n')
  p <- attr(x, 'p')
  if (nargs() > 2) {
    # Case [i, j]
    if (missing(i)) {
      i <- 1:n
    } else if (class(i) == 'logical') {
      i <- which(rep_len(i, n))
    } else if (class(i) == 'character') {
      i <- sapply(i, function (name) {
        which(rownames(x) == name)
      }, USE.NAMES=FALSE)
    }
    if (missing(j)) {
      j <- 1:p
    } else if (class(j) == 'logical') {
      j <- which(rep_len(j, p))
    } else if (class(j) == 'character') {
      j <- sapply(j, function (name) {
        which(colnames(x) == name)
      }, USE.NAMES=FALSE)
    }
    subset <- matrixSubset(x, i, j)
    # Let R handle drop behavior.
    if(drop == TRUE && (nrow(subset) == 1 || ncol(subset) == 1)) {
      subset <- subset[,]
    }
  } else {
    if (missing(i)) {
      # Case []
      i <- 1:n
      j <- 1:p
      subset <- matrixSubset(x, i, j)
    } else {
      # Case [i]
      if (class(i) == 'logical') {
        i <- which(rep_len(i, n * p))
      }
      subset <- vectorSubset(x, i)
    }
  }
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
