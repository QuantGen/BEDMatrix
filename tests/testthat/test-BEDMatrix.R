for (path in c(paste0(extdataPath, "/example"), paste0(extdataPath, "/example.bed"))) {

    test_that("it throws an error if file does not exist", {
        expect_error(BEDMatrix("NOT_FOUND"), "File not found\\.")
    })

    test_that("it throws an error if file is not a BED file", {
        expect_error(BEDMatrix(system.file("extdata", "example.raw", package = "BEDMatrix")), "File is not a binary PED file\\.")
    })

    test_that("it determines n from FAM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(nrow(bed), nrow(CROCHET_EXTRACT_ENV$COMPARE_OBJECT))
        expect_message(BEDMatrix(path = path), "Extracting number of samples and rownames from example\\.fam\\.\\.\\.")
    })

    test_that("it throws an error if FAM file is not found and n is not given", {
        expect_error(BEDMatrix(path = standalonePath), "standalone.fam not found\\. Provide number of samples \\(n\\)\\.")
    })

    test_that("it determines rownames from FAM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(rownames(bed), rownames(CROCHET_EXTRACT_ENV$COMPARE_OBJECT))
        expect_message(BEDMatrix(path = path), "Extracting number of samples and rownames from example\\.fam\\.\\.\\.")
    })

    test_that("it determines p from BIM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(ncol(bed), ncol(CROCHET_EXTRACT_ENV$COMPARE_OBJECT))
        expect_message(BEDMatrix(path = path), "Extracting number of variants and colnames from example\\.bim\\.\\.\\.")
    })

    test_that("it throws an error if BIM file is not found and p is not given", {
        expect_error(BEDMatrix(path = standalonePath), "standalone.fam not found\\. Provide number of samples \\(n\\)\\.")
    })

    test_that("it determines colnames from BIM file", {
        bed <- BEDMatrix(path = path)
        expect_equal(colnames(bed), colnames(CROCHET_EXTRACT_ENV$COMPARE_OBJECT))
        expect_message(BEDMatrix(path = path), "Extracting number of variants and colnames from example\\.bim\\.\\.\\.")
    })

    test_that("it accepts n and p if FAM or BIM file are present", {
        bed <- BEDMatrix(path = path, n = nrow(CROCHET_EXTRACT_ENV$COMPARE_OBJECT), p = ncol(CROCHET_EXTRACT_ENV$COMPARE_OBJECT))
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
