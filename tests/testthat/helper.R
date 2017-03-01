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

extdataPath <- system.file("extdata", package = "BEDMatrix")
# Fix devtools not redirecting system.file calls
if (extdataPath == "") {
    extdataPath <- "../../inst/extdata"
}
standalonePath <- "standalone.bed"

CROCHET_EXTRACT_A <- suppressMessages(BEDMatrix(path = paste0(extdataPath, "/example.bed")))
CROCHET_EXTRACT_B <- parseRaw(paste0(extdataPath, "/example.raw"))

OUT_OF_BOUNDS_INT <- length(CROCHET_EXTRACT_A) + 1
OUT_OF_BOUNDS_CHAR <- "snp1000_U"
