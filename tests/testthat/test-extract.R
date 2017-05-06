extract_tests <- system.file("tests", "testthat", "test-crochet-extract.R", package = "crochet")
if (extract_tests == "") {
    test_that("subsetting", {
        skip("crochet tests have to be installed")
    })
} else {
    source(extract_tests, local = TRUE)
}
