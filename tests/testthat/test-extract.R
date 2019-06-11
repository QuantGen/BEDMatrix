# Source extraction tests
CROCHET_EXTRACT_ENV <- new.env()
CROCHET_EXTRACT_ENV$COMPARE_OBJECT <- raw
CROCHET_EXTRACT_ENV$CUSTOM_OBJECT <- bed
CROCHET_EXTRACT_ENV$OUT_OF_BOUNDS_INT <- length(CROCHET_EXTRACT_ENV$CUSTOM_OBJECT) + 1
CROCHET_EXTRACT_ENV$OUT_OF_BOUNDS_CHAR <- "snp1000_U"
source(system.file("test-suite", "crochet-extract.R", package = "crochet"), local = TRUE)
