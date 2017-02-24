for (path in c(paste0(examplePath, "/example"), paste0(examplePath, "/example.bed"))) {

    test_that("it throws an error if file does not exist", {
        expect_error(BEDMatrix("NOT_FOUND"), "File not found\\.")
    })

    test_that("it throws an error if file is not a BED file", {
        expect_error(BEDMatrix(system.file("extdata", "example.raw", package = "BEDMatrix")), "File is not a binary PED file\\.")
    })

    test_that("it determines n from FAM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(nrow(bed), nrow(TST_B))
        expect_message(BEDMatrix(path = path), "Extracting number of individuals and rownames from FAM file\\.\\.\\.")
    })

    test_that("it throws an error if FAM file is not found and n is not given", {
        expect_error(BEDMatrix(path = standalonePath), "FAM file of same name not found\\. Provide number of individuals \\(n\\)\\.")
    })

    test_that("it determines rownames from FAM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(rownames(bed), rownames(TST_B))
        expect_message(BEDMatrix(path = path), "Extracting number of individuals and rownames from FAM file\\.\\.\\.")
    })

    test_that("it determines p from BIM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(ncol(bed), ncol(TST_B))
        expect_message(BEDMatrix(path = path), "Extracting number of markers and colnames from BIM file\\.\\.\\.")
    })

    test_that("it throws an error if BIM file is not found and p is not given", {
        expect_error(BEDMatrix(path = standalonePath), "FAM file of same name not found\\. Provide number of individuals \\(n\\)\\.")
    })

    test_that("it determines colnames from BIM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(colnames(bed), colnames(TST_B))
        expect_message(BEDMatrix(path = path), "Extracting number of markers and colnames from BIM file\\.\\.\\.")
    })

    test_that("it accepts n and p if FAM or BIM file are present", {
        bed <- BEDMatrix(path = path, n = nrow(TST_B), p = ncol(TST_B))
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
