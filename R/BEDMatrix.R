# Delimiters used in .fam and .bim files
delims <- "[ \t]"

BEDMatrix <- setClass("BEDMatrix", slots = c(
    xptr = "externalptr",
    dims = "integer",
    dnames = "list",
    path = "character"
))

BEDMatrix <- function(path, n = NULL, p = NULL, simple_names = FALSE) {
    path <- path.expand(path)
    if (!file.exists(path)) {
        # Try to add extension (common in PLINK)
        path <- paste0(path, ".bed")
        if (!file.exists(path)) {
            stop("File not found.", call. = FALSE)
        }
    }
    pathSansExt <- tools::file_path_sans_ext(path)
    filesetName <- basename(pathSansExt)
    if (is.null(n)) {
        # Check if .fam file exists
        famPath <- paste0(pathSansExt, ".fam")
        if (!file.exists(famPath)) {
            stop(filesetName, ".fam not found. Provide number of samples (n).", call. = FALSE)
        } else {
            message("Extracting number of samples and rownames from ", filesetName, ".fam...")
            if (requireNamespace("data.table", quietly = TRUE)) {
                if (simple_names) {
                    famColumns <- c(2L)
                } else {
                    famColumns <- c(1L, 2L)
                }
                fam <- data.table::fread(
                    famPath,
                    select = famColumns,
                    colClasses = list(character = famColumns),
                    data.table = FALSE,
                    showProgress = FALSE
                )
                # Determine n
                n <- nrow(fam)
                # Determine rownames
                if (simple_names) {
                    # Use within-family ID only
                    rownames <- fam[, 1L]
                } else {
                    # Concatenate family ID and within-family ID
                    rownames <- paste0(fam[, 1L], "_", fam[, 2L])
                }
            } else {
                fam <- readLines(famPath) # much faster than read.table
                # Determine n
                n <- length(fam)
                # Determine rownames
                if (simple_names) {
                    # Use within-family ID only
                    rownames <- vapply(strsplit(fam, delims), `[`, "", 2L)
                } else {
                    rownames <- vapply(strsplit(fam, delims), function(line) {
                        # Concatenate family ID and within-family ID
                        paste0(line[1L], "_", line[2L])
                    }, "")
                }
            }
        }
    } else {
        n <- as.integer(n)
        rownames <- NULL
    }
    if (is.null(p)) {
        # Check if .bim file exists
        bimPath <- paste0(pathSansExt, ".bim")
        if (!file.exists(bimPath)) {
            stop(filesetName, ".bim not found. Provide number of variants (p).", .call = FALSE)
        } else {
            message("Extracting number of variants and colnames from ", filesetName, ".bim...")
            if (requireNamespace("data.table", quietly = TRUE)) {
                if (simple_names) {
                    bimColumns <- c(2L)
                } else {
                    bimColumns <- c(2L, 5L)
                }
                bim <- data.table::fread(
                    bimPath,
                    select = bimColumns,
                    colClasses = list(character = bimColumns),
                    data.table = FALSE,
                    showProgress = FALSE
                )
                # Determine p
                p <- nrow(bim)
                # Determine colnames
                if (simple_names) {
                    # Use variant name only
                    colnames <- bim[, 1L]
                } else {
                    # Concatenate variant name and A1 allele (like '--recode A'
                    # in PLINK)
                    colnames <- paste0(bim[, 1L], "_", bim[, 2L])
                }
            } else {
                bim <- readLines(bimPath) # much faster than read.table
                # Determine p
                p <- length(bim)
                # Determine colnames
                if (simple_names) {
                    # Use variant name only
                    colnames <- vapply(strsplit(bim, delims), `[`, "", 2L)
                } else {
                    colnames <- vapply(strsplit(bim, delims), function(line) {
                        # Concatenate variant name and A1 allele (like
                        # '--recode A' in PLINK)
                        paste0(line[2L], "_", line[5L])
                    }, "")
                }
            }
        }
    } else {
        p <- as.integer(p)
        colnames <- NULL
    }
    obj <- new(
        "BEDMatrix",
        xptr = .Call(C_BEDMatrix_initialize, path, n, p),
        path = path,
        dims = c(n, p),
        dnames = list(rownames, colnames)
    )
    return(obj)
}

setMethod("show", "BEDMatrix", function(object) {
    dims <- dim(object)
    n <- dims[1L]
    p <- dims[2L]
    cat("BEDMatrix: ", n, " x ", p, " [", slot(object, "path"), "]\n", sep = "")
})

extract_vector <- function(x, i) {
    .Call(C_BEDMatrix_extract_vector, slot(x, "xptr"), i)
}

extract_matrix <- function(x, i, j) {
    subset <- .Call(C_BEDMatrix_extract_matrix, slot(x, "xptr"), i, j)
    # Preserve dimnames
    names <- slot(x, "dnames")
    dimnames(subset) <- list(
        names[[1L]][i],
        names[[2L]][j]
    )
    return(subset)
}

`[.BEDMatrix` <- extract(
    extract_vector = extract_vector,
    extract_matrix = extract_matrix,
    allowDoubles = TRUE
)

dim.BEDMatrix <- function(x) {
    slot(x, "dims")
}

dimnames.BEDMatrix <- function(x) {
    slot(x, "dnames")
}

`dimnames<-.BEDMatrix` <- function(x, value) {
    d <- dim(x)
    v1 <- value[[1L]]
    v2 <- value[[2L]]
    if (!is.list(value) || length(value) != 2L || !(is.null(v1) || length(v1) == d[1L]) || !(is.null(v2) || length(v2) == d[2L])) {
        stop("invalid dimnames", call. = FALSE)
    }
    slot(x, "dnames") <- lapply(value, function(v) {
        if (!is.null(v)) {
            as.character(v)
        }
    })
    return(x)
}

length.BEDMatrix <- function(x) {
    prod(dim(x))
}

str.BEDMatrix <- function(object, ...) {
    print(object)
}

as.matrix.BEDMatrix <- function(x, ...) {
    x[, , drop = FALSE]
}

is.matrix.BEDMatrix <- function(x) {
    TRUE
}
