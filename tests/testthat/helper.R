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

examplePath <- system.file("extdata", package = "BEDMatrix")
# Fix devtools not redirecting system.file calls
if (examplePath == "") {
    examplePath <- "../../inst/extdata"
}
standalonePath <- "standalone.bed"

TST_A <- suppressMessages(BEDMatrix(path = paste0(examplePath, "/example.bed")))
TST_B <- parseRaw(paste0(examplePath, "/example.raw"))

OUT_OF_BOUNDS_INT <- length(TST_A) + 1
OUT_OF_BOUNDS_CHAR <- "snp1000_U"
