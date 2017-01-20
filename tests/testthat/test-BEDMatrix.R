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

test_that("subsetting", {

    expect_equal(bed[], raw[])
    expect_equal(typeof(bed[]), "integer")

    expect_equal(bed[1, ], raw[1, ])
    expect_equal(bed[, 1], raw[, 1])
    expect_equal(bed[1, 1], raw[1, 1])
    expect_equal(bed[1, , drop = FALSE], raw[1, , drop = FALSE])
    expect_equal(bed[, 1, drop = FALSE], raw[, 1, drop = FALSE])
    expect_equal(bed[1, 1, drop = FALSE], raw[1, 1, drop = FALSE])
    expect_equal(typeof(bed[1, ]), "integer")

    expect_equal(bed[1:2, ], raw[1:2, ])
    expect_equal(bed[, 1:2], raw[, 1:2])
    expect_equal(bed[1:2, 1:2], raw[1:2, 1:2])
    expect_equal(bed[1:2, , drop = FALSE], raw[1:2, , drop = FALSE])
    expect_equal(bed[, 1:2, drop = FALSE], raw[, 1:2, drop = FALSE])
    expect_equal(bed[1:2, 1:2, drop = FALSE], raw[1:2, 1:2, drop = FALSE])
    expect_equal(typeof(bed[1:2, ]), "integer")

    expect_equal(bed[2:1, ], raw[2:1, ])
    expect_equal(bed[, 2:1], raw[, 2:1])
    expect_equal(bed[2:1, 2:1], raw[2:1, 2:1])
    expect_equal(bed[2:1, , drop = FALSE], raw[2:1, , drop = FALSE])
    expect_equal(bed[, 2:1, drop = FALSE], raw[, 2:1, drop = FALSE])
    expect_equal(bed[2:1, 2:1, drop = FALSE], raw[2:1, 2:1, drop = FALSE])
    expect_equal(typeof(bed[2:1, ]), "integer")

    expect_equal(bed[c(3, 1), ], raw[c(3, 1), ])
    expect_equal(bed[, c(3, 1)], raw[, c(3, 1)])
    expect_equal(bed[c(3, 1), c(3, 1)], raw[c(3, 1), c(3, 1)])
    expect_equal(bed[c(3, 1), , drop = FALSE], raw[c(3, 1), , drop = FALSE])
    expect_equal(bed[, c(3, 1), drop = FALSE], raw[, c(3, 1), drop = FALSE])
    expect_equal(bed[c(3, 1), c(3, 1), drop = FALSE], raw[c(3, 1), c(3, 1), drop = FALSE])
    expect_equal(typeof(bed[c(3, 1), ]), "integer")

    expect_equal(bed[c(TRUE, FALSE, TRUE), ], raw[c(TRUE, FALSE, TRUE), ])
    expect_equal(bed[, c(TRUE, FALSE, TRUE)], raw[, c(TRUE, FALSE, TRUE)])
    expect_equal(bed[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)], raw[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)])
    expect_equal(bed[c(TRUE, FALSE, TRUE), , drop = FALSE], raw[c(TRUE, FALSE, TRUE), , drop = FALSE])
    expect_equal(bed[, c(TRUE, FALSE, TRUE), drop = FALSE], raw[, c(TRUE, FALSE, TRUE), drop = FALSE])
    expect_equal(bed[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE], raw[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE])
    expect_equal(typeof(bed[c(TRUE, FALSE, TRUE), ]), "integer")

    expect_equal(bed[1], raw[1])
    expect_equal(bed[1:2], raw[1:2])
    expect_equal(bed[2:1], raw[2:1])
    expect_equal(bed[c(3, 1)], raw[c(3, 1)])
    expect_equal(bed[c(TRUE, FALSE, TRUE)], raw[c(TRUE, FALSE, TRUE)])
    expect_equal(bed[raw > 1], raw[raw > 1])
    expect_equal(typeof(bed[1]), "integer")

    expect_equal(bed["per0_per0", ], raw["per0_per0", ])
    expect_equal(bed[, "snp0_A"], raw[, "snp0_A"])
    expect_equal(bed["per0_per0", "snp0_A"], raw["per0_per0", "snp0_A"])
    expect_equal(bed["per0_per0", , drop = FALSE], raw["per0_per0", , drop = FALSE])
    expect_equal(bed[, "snp0_A", drop = FALSE], raw[, "snp0_A", drop = FALSE])
    expect_equal(bed["per0_per0", "snp0_A", drop = FALSE], raw["per0_per0", "snp0_A", drop = FALSE])
    expect_equal(typeof(bed["per0_per0", ]), "integer")

    expect_equal(bed[c("per0_per0", "per1_per1"), ], raw[c("per0_per0", "per1_per1"), ])
    expect_equal(bed[, c("snp0_A", "snp1_C")], raw[, c("snp0_A", "snp1_C")])
    expect_equal(bed[c("per0_per0", "per1_per1"), c("snp0_A", "snp1_C")], raw[c("per0_per0", "per1_per1"), c("snp0_A", "snp1_C")])
    expect_equal(bed[c("per0_per0", "per1_per1"), , drop = FALSE], raw[c("per0_per0", "per1_per1"), , drop = FALSE])
    expect_equal(bed[, c("snp0_A", "snp1_C"), drop = FALSE], raw[, c("snp0_A", "snp1_C"), drop = FALSE])
    expect_equal(bed[c("per0_per0", "per1_per1"), c("snp0_A", "snp1_C"), drop = FALSE], raw[c("per0_per0", "per1_per1"), c("snp0_A", "snp1_C"), drop = FALSE])
    expect_equal(typeof(bed[c("per0_per0", "per1_per1"), ]), "integer")

    expect_equal(bed[c("per1_per1", "per0_per0"), ], raw[c("per1_per1", "per0_per0"), ])
    expect_equal(bed[, c("snp1_C", "snp0_A")], raw[, c("snp1_C", "snp0_A")])
    expect_equal(bed[c("per1_per1", "per0_per0"), c("snp1_C", "snp0_A")], raw[c("per1_per1", "per0_per0"), c("snp1_C", "snp0_A")])
    expect_equal(bed[c("per1_per1", "per0_per0"), , drop = FALSE], raw[c("per1_per1", "per0_per0"), , drop = FALSE])
    expect_equal(bed[, c("snp1_C", "snp0_A"), drop = FALSE], raw[, c("snp1_C", "snp0_A"), drop = FALSE])
    expect_equal(bed[c("per1_per1", "per0_per0"), c("snp1_C", "snp0_A"), drop = FALSE], raw[c("per1_per1", "per0_per0"), c("snp1_C", "snp0_A"), drop = FALSE])
    expect_equal(typeof(bed[c("per1_per1", "per0_per0"), ]), "integer")

    expect_equal(bed[c("per2_per2", "per0_per0"), ], raw[c("per2_per2", "per0_per0"), ])
    expect_equal(bed[, c("snp2_G", "snp0_A")], raw[, c("snp2_G", "snp0_A")])
    expect_equal(bed[c("per2_per2", "per0_per0"), c("snp2_G", "snp0_A")], raw[c("per2_per2", "per0_per0"), c("snp2_G", "snp0_A")])
    expect_equal(bed[c("per2_per2", "per0_per0"), , drop = FALSE], raw[c("per2_per2", "per0_per0"), , drop = FALSE])
    expect_equal(bed[, c("snp2_G", "snp0_A"), drop = FALSE], raw[, c("snp2_G", "snp0_A"), drop = FALSE])
    expect_equal(bed[c("per2_per2", "per0_per0"), c("snp2_G", "snp0_A"), drop = FALSE], raw[c("per2_per2", "per0_per0"), c("snp2_G", "snp0_A"), drop = FALSE])
    expect_equal(typeof(bed[c("per2_per2", "per0_per0"), ]), "integer")

    # Do not modify indexes
    i <- seq_len(nrow(raw) * ncol(raw))
    bed[i]
    expect_equal(i, seq_len(nrow(raw) * ncol(raw)))

    i <- seq_len(nrow(raw))
    j <- seq_len(ncol(raw))
    bed[i, j]
    expect_equal(i, seq_len(nrow(raw)))
    expect_equal(j, seq_len(ncol(raw)))

})

test_that("length", {
    expect_equal(length(bed), length(raw))
})
