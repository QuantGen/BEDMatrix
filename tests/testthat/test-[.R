extract_tests <- system.file("tests", "testthat", "test-[.R", package = "crochet")
if (extract_tests == "") {
    skip("Skip subsetting tests: crochet tests have to be installed")
} else {
    source(extract_tests, local = TRUE)
}
