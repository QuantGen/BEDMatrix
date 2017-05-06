# Delimiters used in PED files
delims <- "[ \t]"

initialize <- function(.Object, path, n = NULL, p = NULL) {
    path <- path.expand(path)
    if (!file.exists(path)) {
        # Try to add extension (common in PLINK)
        path <- paste0(path, ".bed")
        if (!file.exists(path)) {
            stop("File not found.")
        }
    }
    dir <- substr(path, 1L, nchar(path) - 4L)
    if (is.null(n)) {
        # Check if FAM file exists
        famPath <- paste0(dir, ".fam")
        if (!file.exists(famPath)) {
            stop("FAM file of same name not found. Provide number of samples (n).")
        } else {
            message("Extracting number of samples and rownames from FAM file...")
            if (requireNamespace("data.table", quietly = TRUE)) {
                fam <- data.table::fread(famPath, select = c(1L, 2L), data.table = FALSE, showProgress = FALSE)
                # Determine n
                n <- nrow(fam)
                # Determine rownames
                rownames <- paste0(fam[, 1L], "_", fam[, 2L])
            } else {
                fam <- readLines(famPath)
                # Determine n
                n <- length(fam)
                # Determine rownames
                rownames <- sapply(strsplit(fam, delims), function(line) {
                    # Concatenate family ID and subject ID
                    return(paste0(line[1L], "_", line[2L]))
                })
            }
        }
    } else {
        n <- as.integer(n)
        rownames <- NULL
    }
    if (is.null(p)) {
        # Check if BIM file exists
        bimPath <- paste0(dir, ".bim")
        if (!file.exists(bimPath)) {
            stop("BIM file of same name not found. Provide number of variants (p).")
        } else {
            message("Extracting number of variants and colnames from BIM file...")
            if (requireNamespace("data.table", quietly = TRUE)) {
                bim <- data.table::fread(bimPath, select = c(2L, 5L), data.table = FALSE, showProgress = FALSE)
                # Determine p
                p <- nrow(bim)
                # Determine colnames
                colnames <- paste0(bim[, 1L], "_", bim[, 2L])
            } else {
                bim <- readLines(bimPath)
                # Determine p
                p <- length(bim)
                # Determine colnames
                colnames <- sapply(strsplit(bim, delims), function(line) {
                    # Concatenate SNP name and minor allele (like --recodeA)
                    return(paste0(line[2L], "_", line[5L]))
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

extract_vector <- function(x, i) {
    .Call("BEDMatrix__extract_vector", x@xptr, i)
}

extract_matrix <- function(x, i, j) {
    subset <- .Call("BEDMatrix__extract_matrix", x@xptr, i, j)
    # Preserve dimnames
    names <- x@dnames
    dimnames(subset) <- list(
        names[[1L]][i],
        names[[2L]][j]
    )
    return(subset)
}

show <- function(object) {
    dims <- dim(object)
    n <- dims[1L]
    p <- dims[2L]
    cat("BEDMatrix: ", n, " x ", p, " [", object@path, "]\n", sep = "")
}

#' A Class to Extract Genotypes from a PLINK .bed File.
#'
#' `BEDMatrix` is a class that maps a [PLINK .bed](https://www.cog-genomics.org/plink2/formats#bed)
#' file into memory and behaves similarly to a regular `matrix` by implementing
#' key methods such as `[`, `dim`, and `dimnames`. Subsets are extracted
#' directly and on-demand from the .bed file without loading the entire file
#' into memory.
#'
#' The subsets extracted from a `BEDMatrix` object are coded similarly to
#' [.raw](https://www.cog-genomics.org/plink2/formats#raw) files (generated
#' with the `--recodeA` argument in
#' [PLINK](https://www.cog-genomics.org/plink2/)): `0` indicates homozygous
#' major allele, `1` indicates heterozygous, and `2` indicates homozygous minor
#' allele.
#'
#' Internally, this class is an S4 class with the following slots that should
#' not be relied upon in actual code: `xptr`, `dims`, `dnames`, and `path`. The
#' .bed file is mapped into memory using the [Rcpp][Rcpp::Rcpp-package] package
#' and the `Boost.Interprocess` library provided by the [BH][BH::BH-package]
#' package.
#'
#' @section Methods:
#' - `[`
#' - `dim`
#' - `dimnames`
#' - `dimnames<-`
#' - `as.matrix`
#' - `is.matrix`
#' - `length`
#' - `str`
#' - `show`
#' - `initialize`
#'
#' @slot xptr An external pointer to the underlying [Rcpp][Rcpp::Rcpp-package]
#' code.
#' @slot dims An integer vector specifying the number of samples and variants
#' as determined by the the accompanying
#' [.fam](https://www.cog-genomics.org/plink2/formats#fam) and
#' [.bim](https://www.cog-genomics.org/plink2/formats#bim) files or by the `n`
#' and `p` parameters of the [constructor
#' function][initialize,BEDMatrix-method()].
#' @slot dnames A list describing the row names and column names of the object
#' as determined by the accompanying
#' [.fam](https://www.cog-genomics.org/plink2/formats#fam) and
#' [.bim](https://www.cog-genomics.org/plink2/formats#bim) files, or `NULL` if
#' the `n` and `p` parameters of the [constructor
#' function][initialize,BEDMatrix-method()] were provided.
#' @slot path A character string containing the path to the .bed file.
#' @seealso [initialize()][initialize,BEDMatrix-method()] to create a
#' `BEDMatrix` object from a .bed file, [BEDMatrix-package] to learn more about
#' .bed files, [LinkedMatrix][LinkedMatrix::LinkedMatrix-package] to link
#' several `BEDMatrix` objects together.
#' @example man/examples/BEDMatrix.R
#' @aliases BEDMatrix-class
#' @export BEDMatrix
#' @exportClass BEDMatrix
BEDMatrix <- setClass("BEDMatrix", slots = c(xptr = "externalptr", dims = "integer", dnames = "list", path = "character"))

#' Create a BEDMatrix Object from a PLINK .bed File.
#'
#' This function constructs a new [BEDMatrix-class] object by mapping the
#' specified [PLINK .bed](https://www.cog-genomics.org/plink2/formats#bed) file
#' into memory.
#'
#' [.bed](https://www.cog-genomics.org/plink2/formats#bed) files must be
#' accompanied by [.fam](https://www.cog-genomics.org/plink2/formats#fam) and
#' [.bim](https://www.cog-genomics.org/plink2/formats#bim) files: .fam files
#' contain sample information, and .bim files contain variant information. If
#' the name of the .bed file is *plink*.bed then the names of the .fam and .bim
#' files have to be *plink*.fam and *plink*.bim, respectively. The .fam and
#' .bim files are used to extract the number and names of samples and variants.
#'
#' For very large .bed files, reading the .fam and .bim files can take a long
#' time. If `n` and `p` are provided, these files are not read and `dimnames`
#' have to be provided manually.
#'
#' Currently, only the variant-major mode of .bed files is supported.
#' [PLINK2](https://www.cog-genomics.org/plink2/) "dropped" support for the
#' sample-major mode by automatically converting files in this format to the
#' variant-major mode. Therefore, it is recommended to run files in
#' sample-major mode through PLINK2 first.
#'
#' @param .Object Internal, used by [methods::initialize()] generic.
#' @param path Path to the
#' [.bed](https://www.cog-genomics.org/plink2/formats#bed) file (with or
#' without extension).
#' @param n The number of samples. If `NULL` (the default), this number will be
#' determined from the accompanying
#' [.fam](https://www.cog-genomics.org/plink2/formats#fam) file (of the same
#' name as the [.bed](https://www.cog-genomics.org/plink2/formats#bed) file).
#' If a positive integer, the .fam file is not read and `rownames` will be set
#' to `NULL` and have to be provided manually.
#' @param p The number of variants. If `NULL` (the default) the number of
#' variants will be determined from the accompanying
#' [.bim](https://www.cog-genomics.org/plink2/formats#bim) file (of the same
#' name as the [.bed](https://www.cog-genomics.org/plink2/formats#bed) file).
#' If a positive integer, the .bim file is not read and `colnames` will be set
#' to `NULL` and have to be provided manually.
#' @return A [BEDMatrix-class] object.
#' @example man/examples/initialize.R
#' @seealso [BEDMatrix-package] to learn more about .bed files.
#' @export
setMethod("initialize", signature(.Object = "BEDMatrix"), initialize)

#' Show a BEDMatrix Object.
#'
#' Display the object, by printing, plotting or whatever suits its class.
#'
#' @param object A [BEDMatrix-class] object.
#' @export
setMethod("show", signature(object = "BEDMatrix"), show)

#' @export
`[.BEDMatrix` <- crochet::extract(extract_vector = extract_vector, extract_matrix = extract_matrix)

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
    v1 <- value[[1L]]
    v2 <- value[[2L]]
    if (!is.list(value) || length(value) != 2L || !(is.null(v1) || length(v1) == d[1L]) || !(is.null(v2) || length(v2) == d[2L])) {
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
as.matrix.BEDMatrix <- function(x, ...) {
    x[, , drop = FALSE]
}

#' @export
is.matrix.BEDMatrix <- function(x) {
    TRUE
}
