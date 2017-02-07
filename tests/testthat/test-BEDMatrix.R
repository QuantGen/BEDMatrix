context("BEDMatrix")

parseRaw <- function(path) {
    lines <- strsplit(readLines(path), " ")
    header <- lines[[1]]
    data <- matrix(data = unlist(lines[2:length(lines)]), nrow = 50, ncol = 1006, byrow = TRUE)
    pheno <- data[, 1:6]
    geno <- data[, 7:ncol(data)]
    suppressWarnings(mode(geno) <- "integer")
    rownames(geno) <- paste0(pheno[, 1], "_", pheno[, 2])
    colnames(geno) <- header[7:length(header)]
    return(geno)
}

# Prepare dummy data
raw <- parseRaw(system.file("extdata", "example.raw", package = "BEDMatrix"))
examplePath <- system.file("extdata", "example.bed", package = "BEDMatrix")
standalonePath <- "standalone.bed"

for (path in c(examplePath, sub(".bed", "", examplePath))) {

    test_that("it throws an error if file does not exist", {
        expect_error(BEDMatrix("NOT_FOUND"), "File not found\\.")
    })

    test_that("it throws an error if file is not a BED file", {
        expect_error(BEDMatrix(system.file("extdata", "example.raw", package = "BEDMatrix")), "File is not a binary PED file\\.")
    })

    test_that("it determines n from FAM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(nrow(bed), nrow(raw))
        expect_message(BEDMatrix(path = path), "Extracting number of individuals and rownames from FAM file\\.\\.\\.")
    })

    test_that("it throws an error if FAM file is not found and n is not given", {
        expect_error(BEDMatrix(path = standalonePath), "FAM file of same name not found\\. Provide number of individuals \\(n\\)\\.")
    })

    test_that("it determines rownames from FAM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(rownames(bed), rownames(raw))
        expect_message(BEDMatrix(path = path), "Extracting number of individuals and rownames from FAM file\\.\\.\\.")
    })

    test_that("it determines p from BIM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(ncol(bed), ncol(raw))
        expect_message(BEDMatrix(path = path), "Extracting number of markers and colnames from BIM file\\.\\.\\.")
    })

    test_that("it throws an error if BIM file is not found and p is not given", {
        expect_error(BEDMatrix(path = standalonePath), "FAM file of same name not found\\. Provide number of individuals \\(n\\)\\.")
    })

    test_that("it determines colnames from BIM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(colnames(bed), colnames(raw))
        expect_message(BEDMatrix(path = path), "Extracting number of markers and colnames from BIM file\\.\\.\\.")
    })

    test_that("it accepts n and p if FAM or BIM file are present", {
        bed <- BEDMatrix(path = path, n = nrow(raw), p = ncol(raw))
        expect_equal(dimnames(bed), list(NULL, NULL))
    })

    test_that("it accepts n and p if FAM or BIM file is not found", {
        bed <- BEDMatrix(path = standalonePath, n = 3, p = 6)
        expect_equal(dimnames(bed), list(NULL, NULL))
    })

    test_that("it throws an error if dimensions are wrong", {
        expect_error(BEDMatrix(path = path, n = 10, p = 5), "n or p does not match the dimensions of the file\\.")
    })

}

# Prepare dummy BED matrix
suppressMessages(bed <- BEDMatrix(path = examplePath))

test_that("length", {
    expect_equal(length(bed), length(raw))
})
