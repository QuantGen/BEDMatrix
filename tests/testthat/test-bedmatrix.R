context("BEDMatrix")

parseRaw <- function(path) {
    file <- file(path)
    lines <- readLines(file)
    out <- t(sapply(strsplit(lines[2:length(lines)], " "), function(line) {
        line <- line[7:length(line)]
        line <- gsub("NA", NA, line)
        as.integer(line)
    }))
    rownames(out) <- paste0("id_", 1:nrow(out))
    colnames(out) <- paste0("mrk_", 1:ncol(out))
    close(file)
    return(out)
}

# Prepare dummy data
raw <- parseRaw(system.file("extdata", "example.raw", package = "BEDMatrix"))

test_that("it throws an error if file does not exist", {
    expect_error(bed("NOT_FOUND"))
})

test_that("it determines n from FAM file", {
    bed <- BEDMatrix(path = system.file("extdata", "example.bed", package = "BEDMatrix"))
    expect_equal(nrow(bed), 6)
})

test_that("subsetting", {

    bed <- BEDMatrix(path = system.file("extdata", "example.bed", package = "BEDMatrix"), n = 6, p = 3)
    rownames(bed) <- paste0("id_", 1:nrow(bed))
    colnames(bed) <- paste0("mrk_", 1:ncol(bed))

    expect_true(all.equal(bed[], raw[]))

    expect_true(all.equal(bed[1, ], raw[1, ]))
    expect_true(all.equal(bed[, 1], raw[, 1]))
    expect_true(all.equal(bed[1, 1], raw[1, 1]))
    expect_true(all.equal(bed[1, , drop = FALSE], raw[1, , drop = FALSE]))
    expect_true(all.equal(bed[, 1, drop = FALSE], raw[, 1, drop = FALSE]))
    expect_true(all.equal(bed[1, 1, drop = FALSE], raw[1, 1, drop = FALSE]))

    expect_true(all.equal(bed[1:2, ], raw[1:2, ]))
    expect_true(all.equal(bed[, 1:2], raw[, 1:2]))
    expect_true(all.equal(bed[1:2, 1:2], raw[1:2, 1:2]))
    expect_true(all.equal(bed[1:2, , drop = FALSE], raw[1:2, , drop = FALSE]))
    expect_true(all.equal(bed[, 1:2, drop = FALSE], raw[, 1:2, drop = FALSE]))
    expect_true(all.equal(bed[1:2, 1:2, drop = FALSE], raw[1:2, 1:2, drop = FALSE]))

    expect_true(all.equal(bed[2:1, ], raw[2:1, ]))
    expect_true(all.equal(bed[, 2:1], raw[, 2:1]))
    expect_true(all.equal(bed[2:1, 2:1], raw[2:1, 2:1]))
    expect_true(all.equal(bed[2:1, , drop = FALSE], raw[2:1, , drop = FALSE]))
    expect_true(all.equal(bed[, 2:1, drop = FALSE], raw[, 2:1, drop = FALSE]))
    expect_true(all.equal(bed[2:1, 2:1, drop = FALSE], raw[2:1, 2:1, drop = FALSE]))

    expect_true(all.equal(bed[c(3, 1), ], raw[c(3, 1), ]))
    expect_true(all.equal(bed[, c(3, 1)], raw[, c(3, 1)]))
    expect_true(all.equal(bed[c(3, 1), c(3, 1)], raw[c(3, 1), c(3, 1)]))
    expect_true(all.equal(bed[c(3, 1), , drop = FALSE], raw[c(3, 1), , drop = FALSE]))
    expect_true(all.equal(bed[, c(3, 1), drop = FALSE], raw[, c(3, 1), drop = FALSE]))
    expect_true(all.equal(bed[c(3, 1), c(3, 1), drop = FALSE], raw[c(3, 1), c(3, 1), drop = FALSE]))

    expect_true(all.equal(bed[c(TRUE, FALSE, TRUE), ], raw[c(TRUE, FALSE, TRUE), ]))
    expect_true(all.equal(bed[, c(TRUE, FALSE, TRUE)], raw[, c(TRUE, FALSE, TRUE)]))
    expect_true(all.equal(bed[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)], raw[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)]))
    expect_true(all.equal(bed[c(TRUE, FALSE, TRUE), , drop = FALSE], raw[c(TRUE, FALSE, TRUE), , drop = FALSE]))
    expect_true(all.equal(bed[, c(TRUE, FALSE, TRUE), drop = FALSE], raw[, c(TRUE, FALSE, TRUE), drop = FALSE]))
    expect_true(all.equal(bed[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE], raw[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE]))

    expect_true(all.equal(bed[1], raw[1]))
    expect_true(all.equal(bed[1:2], raw[1:2]))
    expect_true(all.equal(bed[2:1], raw[2:1]))
    expect_true(all.equal(bed[c(3, 1)], raw[c(3, 1)]))
    expect_true(all.equal(bed[c(TRUE, FALSE, TRUE)], raw[c(TRUE, FALSE, TRUE)]))
    expect_true(all.equal(bed[raw > 1], raw[raw > 1]))

    expect_true(all.equal(bed["id_1", ], raw["id_1", ]))
    expect_true(all.equal(bed[, "mrk_1"], raw[, "mrk_1"]))
    expect_true(all.equal(bed["id_1", "mrk_1"], raw["id_1", "mrk_1"]))
    expect_true(all.equal(bed["id_1", , drop = FALSE], raw["id_1", , drop = FALSE]))
    expect_true(all.equal(bed[, "mrk_1", drop = FALSE], raw[, "mrk_1", drop = FALSE]))
    expect_true(all.equal(bed["id_1", "mrk_1", drop = FALSE], raw["id_1", "mrk_1", drop = FALSE]))

    expect_true(all.equal(bed[c("id_1", "id_2"), ], raw[c("id_1", "id_2"), ]))
    expect_true(all.equal(bed[, c("mrk_1", "mrk_2")], raw[, c("mrk_1", "mrk_2")]))
    expect_true(all.equal(bed[c("id_1", "id_2"), c("mrk_1", "mrk_2")], raw[c("id_1", "id_2"), c("mrk_1", "mrk_2")]))
    expect_true(all.equal(bed[c("id_1", "id_2"), , drop = FALSE], raw[c("id_1", "id_2"), , drop = FALSE]))
    expect_true(all.equal(bed[, c("mrk_1", "mrk_2"), drop = FALSE], raw[, c("mrk_1", "mrk_2"), drop = FALSE]))
    expect_true(all.equal(bed[c("id_1", "id_2"), c("mrk_1", "mrk_2"), drop = FALSE], raw[c("id_1", "id_2"), c("mrk_1", "mrk_2"), drop = FALSE]))

    expect_true(all.equal(bed[c("id_2", "id_1"), ], raw[c("id_2", "id_1"), ]))
    expect_true(all.equal(bed[, c("mrk_2", "mrk_1")], raw[, c("mrk_2", "mrk_1")]))
    expect_true(all.equal(bed[c("id_2", "id_1"), c("mrk_2", "mrk_1")], raw[c("id_2", "id_1"), c("mrk_2", "mrk_1")]))
    expect_true(all.equal(bed[c("id_2", "id_1"), , drop = FALSE], raw[c("id_2", "id_1"), , drop = FALSE]))
    expect_true(all.equal(bed[, c("mrk_2", "mrk_1"), drop = FALSE], raw[, c("mrk_2", "mrk_1"), drop = FALSE]))
    expect_true(all.equal(bed[c("id_2", "id_1"), c("mrk_2", "mrk_1"), drop = FALSE], raw[c("id_2", "id_1"), c("mrk_2", "mrk_1"), drop = FALSE]))

    expect_true(all.equal(bed[c("id_3", "id_1"), ], raw[c("id_3", "id_1"), ]))
    expect_true(all.equal(bed[, c("mrk_3", "mrk_1")], raw[, c("mrk_3", "mrk_1")]))
    expect_true(all.equal(bed[c("id_3", "id_1"), c("mrk_3", "mrk_1")], raw[c("id_3", "id_1"), c("mrk_3", "mrk_1")]))
    expect_true(all.equal(bed[c("id_3", "id_1"), , drop = FALSE], raw[c("id_3", "id_1"), , drop = FALSE]))
    expect_true(all.equal(bed[, c("mrk_3", "mrk_1"), drop = FALSE], raw[, c("mrk_3", "mrk_1"), drop = FALSE]))
    expect_true(all.equal(bed[c("id_3", "id_1"), c("mrk_3", "mrk_1"), drop = FALSE], raw[c("id_3", "id_1"), c("mrk_3", "mrk_1"), drop = FALSE]))

})
