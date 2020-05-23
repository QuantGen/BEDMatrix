source("setup.R")

# it throws an error if file does not exist
expect_error(BEDMatrix("NOT_FOUND"), "File not found\\.")

# it throws an error if file is not a BED file
expect_error(BEDMatrix("test-BEDMatrix.R"))

# test both prefix and .bed paths
for (path in c(paste0(extdataPath, "/example"), paste0(extdataPath, "/example.bed"))) {

    # it determines n from FAM file
    bed <- suppressMessages(BEDMatrix(path = path))
    expect_equal(nrow(bed), nrow(raw))
    expect_message(BEDMatrix(path = path), "Extracting number of samples and rownames from example\\.fam\\.\\.\\.")

    # it throws an error if FAM file is not found and n is not given
    expect_error(BEDMatrix(path = standalonePath), "standalone.fam not found\\. Provide number of samples \\(n\\)\\.")

    # it determines rownames from FAM file
    bed <- suppressMessages(BEDMatrix(path = path))
    expect_equal(rownames(bed), rownames(raw))
    expect_message(BEDMatrix(path = path), "Extracting number of samples and rownames from example\\.fam\\.\\.\\.")

    # it determines p from BIM file
    bed <- suppressMessages(BEDMatrix(path = path))
    expect_equal(ncol(bed), ncol(raw))
    expect_message(BEDMatrix(path = path), "Extracting number of variants and colnames from example\\.bim\\.\\.\\.")

    # it throws an error if BIM file is not found and p is not given
    expect_error(BEDMatrix(path = standalonePath), "standalone.fam not found\\. Provide number of samples \\(n\\)\\.")

    # it determines colnames from BIM file
    bed <- suppressMessages(BEDMatrix(path = path))
    expect_equal(colnames(bed), colnames(raw))
    expect_message(BEDMatrix(path = path), "Extracting number of variants and colnames from example\\.bim\\.\\.\\.")

    # it accepts n and p if FAM or BIM file are present
    bed <- BEDMatrix(path = path, n = nrow(raw), p = ncol(raw))
    expect_equal(dimnames(bed), list(NULL, NULL))

    # it accepts n and p if FAM or BIM file is not found
    bed <- BEDMatrix(path = standalonePath, n = 3, p = 6)
    expect_equal(dimnames(bed), list(NULL, NULL))

    # it throws an error if dimensions are wrong
    expect_error(BEDMatrix(path = path, n = 10, p = 5), "n or p does not match the dimensions of the file\\.")

}
