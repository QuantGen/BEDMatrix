subsette_tests <- system.file("tests", "testthat", "test-[.R", package = "subsette")
if (subsette_tests == "") {
    skip("Skip subsetting tests: subsette tests have to be installed")
} else {
    source(subsette_tests, local = TRUE)
}
