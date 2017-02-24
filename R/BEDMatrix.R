#' @useDynLib BEDMatrix BEDMatrix__new
initialize <- function(.Object, path, n = NULL, p = NULL) {
    path <- path.expand(path)
    if (!file.exists(path)) {
        # Try to add extension (common in PLINK)
        path <- paste0(path, ".bed")
        if (!file.exists(path)) {
            stop("File not found.")
        }
    }
    dir <- substr(path, 1, nchar(path) - 4)
    if (is.null(n)) {
        # Check if FAM file exists
        if (!file.exists(paste0(dir, ".fam"))) {
            stop("FAM file of same name not found. Provide number of individuals (n).")
        } else {
            message("Extracting number of individuals and rownames from FAM file...")
            famPath <- paste0(dir, ".fam")
            if (requireNamespace("data.table", quietly = TRUE)) {
                fam <- data.table::fread(famPath, select = c(1, 2), data.table = FALSE, showProgress = FALSE)
                # Determine n
                n <- nrow(fam)
                # Determine rownames
                rownames <- paste0(fam[, 1], "_", fam[, 2])
            } else {
                fam <- readLines(famPath)
                # Determine n
                n <- length(fam)
                # Determine rownames
                rownames <- sapply(strsplit(fam, delims), function(line) {
                    # Concatenate family ID and subject ID
                    return(paste0(line[1], "_", line[2]))
                })
            }
        }
    } else {
        n <- as.integer(n)
        rownames <- NULL
    }
    if (is.null(p)) {
        # Check if BIM file exists
        if (!file.exists(paste0(dir, ".bim"))) {
            stop("BIM file of same name not found. Provide number of markers (p).")
        } else {
            message("Extracting number of markers and colnames from BIM file...")
            bimPath <- paste0(dir, ".bim")
            if (requireNamespace("data.table", quietly = TRUE)) {
                bim <- data.table::fread(bimPath, select = c(2, 5), data.table = FALSE, showProgress = FALSE)
                # Determine p
                p <- nrow(bim)
                # Determine colnames
                colnames <- paste0(bim[, 1], "_", bim[, 2])
            } else {
                bim <- readLines(bimPath)
                # Determine p
                p <- length(bim)
                # Determine colnames
                colnames <- sapply(strsplit(bim, delims), function(line) {
                    # Concatenate SNP name and minor allele (like --recodeA)
                    return(paste0(line[2], "_", line[5]))
                })
            }
        }
    } else {
        p <- as.integer(p)
        colnames <- NULL
    }
    # Create Rcpp object
    .Object@xptr <- .Call("BEDMatrix__new", path, n, p)
    .Object@path <- path
    .Object@dims <- c(n, p)
    .Object@dnames <- list(rownames, colnames)
    return(.Object)
}

#' @useDynLib BEDMatrix BEDMatrix__subset_vector
subset_vector <- function(x, i) {
    .Call("BEDMatrix__subset_vector", x@xptr, i)
}

#' @useDynLib BEDMatrix BEDMatrix__subset_matrix
subset_matrix <- function(x, i, j) {
    subset <- .Call("BEDMatrix__subset_matrix", x@xptr, i, j)
    # Preserve dimnames
    names <- x@dnames
    dimnames(subset) <- list(
        names[[1]][i],
        names[[2]][j]
    )
    return(subset)
}

show <- function(object) {
    dims <- dim(object)
    n <- dims[1]
    p <- dims[2]
    cat("BEDMatrix: ", n, " x ", p, " [", object@path, "]\n", sep = "")
}

# Delimiters used in PED files.
delims <- "[ \t]"

#' Creates a Matrix Wrapper Around Binary PED Files.
#'
#' [BEDMatrix-class] is an S3 class that behaves similarly to a regular
#' `matrix` by implementing key methods such as `[`, `dim`, and `dimnames`.
#' Subsets are extracted directly and on-demand from the binary PED file
#' without loading the entire file into memory through memory mapping. The
#' subsets are coded similarly to RAW files generated with the `--recodeA`
#' argument in PLINK: `0` indicates homozygous major allele, `1` indicates
#' heterozygous, and `2` indicates homozygous minor allele.
#'
#' A [BEDMatrix-class] instance can be created by providing the path to the BED
#' file (with or without extension) as `path`, the number of individuals as
#' `n`, and the number of markers as `p`. If a FAM file (which corresponds to
#' the first six columns of a PED file) of the same name and in the same
#' directory as the BED file exists, it is optional to provide `n` and the
#' number of individuals as well as the rownames of the [BEDMatrix-class]
#' object will be detected automatically.  The rownames will be generated based
#' on the IID and FID of each individual, concatenated by `_`.  If a BIM file
#' (which corresponds to the MAP file that accompanies a PED file) of the same
#' name and in the same directory as the BED file exists, it is optional to
#' provide `p` and the number of markers as well as the colnames of the
#' [BEDMatrix-class] object will be detected automatically. The colnames will
#' be generated based on the SNP name and the minor allele, concatenated by `_`
#' (similar to the colnames in a RAW file). For very large BED file it is
#' advised to provide `n` and `p` manually to speed up object creation. In that
#' case `rownames` and `colnames` will be set to `NULL` and have to be
#' specified manually.
#'
#' A BED file can be created from a PED file with
#' [PLINK](https://www.cog-genomics.org/plink2) using `plink --file myfile
#' --make-bed`. BED files are storage and query efficient, and can be
#' transformed back into the original PED file with PLINK using `plink --bfile
#' myfile --recode`.
#'
#' Internally, [BEDMatrix-class] is an S4 class with the following slots that
#' should not be relied upon in actual code: `path`, `dims`, `dnames`, and
#' `xptr`. `path` stores the path to the BED file. `dims` and `dnames` contain
#' the dimensions and dimnames of the [BEDMatrix-class] object. `xptr` is an
#' external pointer that points to the underlying `Rcpp` functions that
#' memory-map the BED file via `Boost.Interprocess` of the `BH` package.
#'
#' @param path Path to the binary PED file, with or without extension.
#' @param n The number of individuals. Optional if FAM file of same name as BED
#' file exists. If provided, `rownames` will be set to `NULL` and have to be
#' provided manually.
#' @param p The number of markers. Optional if BIM file of same name as BED
#' file exists. If provided, `colnames` will be set to `NULL` and have to be
#' provided manually.
#' @return A `BEDMatrix` object.
#' @export BEDMatrix
#' @exportClass BEDMatrix
#' @example man/examples/BEDMatrix.R
BEDMatrix <- setClass("BEDMatrix", slots = list(xptr = "externalptr", dims = "integer", dnames = "list", path = "character"))

#' @export
setMethod("initialize", signature(.Object = "BEDMatrix"), initialize)

#' @export
setMethod("show", signature(object = "BEDMatrix"), show)

#' @export
`[.BEDMatrix` <- crochet::extract(subset_vector = subset_vector, subset_matrix = subset_matrix)

#' @export
dim.BEDMatrix <- function(x) {
    x@dims
}

#' @export
dimnames.BEDMatrix <- function(x) {
    x@dnames
}

#' @export
`dimnames<-.BEDMatrix` <- function(x, value) {
    d <- dim(x)
    v1 <- value[[1]]
    v2 <- value[[2]]
    if (!is.list(value) || length(value) != 2 || !(is.null(v1) || length(v1) == d[1]) || !(is.null(v2) || length(v2) == d[2])) {
        stop("invalid dimnames")
    }
    x@dnames <- lapply(value, function(v) {
        if (!is.null(v)) {
            as.character(v)
        }
    })
    return(x)
}

#' @export
length.BEDMatrix <- function(x) {
    prod(dim(x))
}

#' @export
str.BEDMatrix <- function(object, ...) {
    print(object)
}

#' @export
`[.BEDMatrix` <- crochet::extract(subset_vector = subset_vector, subset_matrix = subset_matrix)

#' @export
as.matrix.BEDMatrix <- function(x, ...) {
    x[, , drop = FALSE]
}

#' @export
is.matrix.BEDMatrix <- function(x) {
    TRUE
}
