source("setup.R")

# Source extraction tests
extractionTests <- new.env()
extractionTests$COMPARE_OBJECT <- raw
extractionTests$CUSTOM_OBJECT <- bed
extractionTests$OUT_OF_BOUNDS_INT <- length(extractionTests$CUSTOM_OBJECT) + 1
extractionTests$OUT_OF_BOUNDS_CHAR <- "snp1000_U"
source(
    file = system.file("test-suite", "crochet-extract.R", package = "crochet"),
    local = extractionTests
)
