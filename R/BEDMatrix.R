# Delimiters used in PED files.
delims <- "[ \t]"

#' Creates a matrix wrapper around a binary PED file.
#'
#' The binary PED file will not be loaded into memory. Rather, subsets are
#' extracted 'on the fly'.
#'
#' If a FAM or MAP file with the same name as the binary PED file (minus the
#' extension) exists within the same folder, \code{n} and \code{p} as well as
#' rownames and colnames will be automatically determined from those files.
#' To skip this, provide \code{n} and \code{p} manually.
#'
#' @param path Path to the binary PED file.
#' @param n The number of individuals. Optional if FAM file of same name as BED
#'   file exists.
#' @param p The number of markers. Optional if MAP file of same name as BED file
#'   exists.
#' @export
BEDMatrix <- function(path, n = NULL, p = NULL) {
    path <- path.expand(path)
    if (!file.exists(path)) {
        stop("File not found.")
    }
    dir <- substr(path, 1, nchar(path) - 4)
    if (is.null(n)) {
        # Check if FAM file exists.
        if (!file.exists(paste0(dir, ".fam"))) {
            stop("FAM file of same name not found. Provide number of individuals (n).")
        } else {
            message("Extracting number of individuals and rownames from FAM file...")
            fam <- readLines(paste0(dir, ".fam"))
            # Determine n.
            n <- length(fam)
            # Determine rownames.
            rownames <- sapply(strsplit(fam, delims), function(line) {
                # Concatenate family ID and subject ID
                return(paste0(line[1], "_", line[2]))
            })
        }
    } else {
        n <- as.integer(n)
        rownames <- NULL
    }
    if (is.null(p)) {
        # Check if BIM file exists.
        if (!file.exists(paste0(dir, ".bim"))) {
            stop("BIM file of same name not found. Provide number of markers (p).")
        } else {
            message("Extracting number of markers and colnames from BIM file...")
            bim <- readLines(paste0(dir, ".bim"))
            # Determine p.
            p <- length(bim)
            # Determine colnames.
            colnames <- sapply(strsplit(bim, delims), function(line) {
                # Concatenate SNP name and minor allele (like --recodeA)
                return(paste0(line[2], "_", line[5]))
            })
        }
    } else {
        p <- as.integer(p)
        colnames <- NULL
    }
    # Create Rcpp object
    rcpp_obj <- new(BEDMatrix_, path, n, p)
    # Wrap object in S3 class
    s3_obj <- list()
    class(s3_obj) <- "BEDMatrix"
    attr(s3_obj, "_instance") <- rcpp_obj
    attr(s3_obj, "dnames") <- list(rownames, colnames)
    return(s3_obj)
}

#' @export
print.BEDMatrix <- function(x, ...) {
    dims <- dim(x)
    n <- dims[1]
    p <- dims[2]
    cat(paste(n, "x", p, "BEDMatrix"), "\n")
}

#' @export
`[.BEDMatrix` <- function(x, i, j, drop = TRUE) {
    rcpp_obj <- attr(x, "_instance")
    dims <- dim(x)
    n <- dims[1]
    p <- dims[2]
    if (nargs() > 2) {
        # Case [i, j]
        if (missing(i)) {
            i <- 1:n
        } else if (class(i) == "logical") {
            i <- which(rep_len(i, n))
        } else if (class(i) == "character") {
            i <- sapply(i, function(name) {
                which(rownames(x) == name)
            }, USE.NAMES = FALSE)
        }
        if (missing(j)) {
            j <- 1:p
        } else if (class(j) == "logical") {
            j <- which(rep_len(j, p))
        } else if (class(j) == "character") {
            j <- sapply(j, function(name) {
                which(colnames(x) == name)
            }, USE.NAMES = FALSE)
        }
        subset <- rcpp_obj$matrixSubset(x, i, j)
        # Let R handle drop behavior.
        if (drop == TRUE && (nrow(subset) == 1 || ncol(subset) == 1)) {
            subset <- subset[, ]
        }
    } else {
        if (missing(i)) {
            # Case []
            i <- 1:n
            j <- 1:p
            subset <- rcpp_obj$matrixSubset(x, i, j)
        } else {
            # Case [i]
            if (class(i) == "matrix") {
                i <- as.vector(i)
                if (class(i) == "logical") {
                  i <- which(rep_len(i, n * p))
                  # matrix treats NAs as TRUE
                  i <- sort(c(i, which(is.na(x[]))))
                }
            } else {
                if (class(i) == "logical") {
                  i <- which(rep_len(i, n * p))
                }
            }
            subset <- rcpp_obj$vectorSubset(x, i)
        }
    }
    return(subset)
}

#' @export
dim.BEDMatrix <- function(x) {
    rcpp_obj <- attr(x, "_instance")
    n <- rcpp_obj$n
    p <- rcpp_obj$p
    return(c(n, p))
}

#' @export
dimnames.BEDMatrix <- function(x) {
    dnames <- attr(x, "dnames")
    return(dnames)
}

#' @export
`dimnames<-.BEDMatrix` <- function(x, value) {
    d <- dim(x)
    v1 <- value[[1]]
    v2 <- value[[2]]
    if (!is.list(value) || length(value) != 2 || !(is.null(v1) || length(v1) == d[1]) || !(is.null(v2) || length(v2) == d[2])) {
        stop("invalid dimnames")
    }
    attr(x, "dnames") <- lapply(value, function(v) {
        if (!is.null(v)) {
            as.character(v)
        }
    })
    return(x)
}
