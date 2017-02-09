crochet_tests <- system.file("tests", "testthat", "test-[.R", package = "crochet")
if (crochet_tests == "") {
    skip("Skip subsetting tests: crochet tests have to be installed")
} else {
    source(crochet_tests, local = TRUE)
}
