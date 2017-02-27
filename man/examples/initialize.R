# Get the path to the example .bed file
path <- system.file("extdata", "example.bed", package = "BEDMatrix")

# Create a BEDMatrix object the example .bed file
m1 <- BEDMatrix(path)

# Create a BEDMatrix object the example .bed file without loading the .fam and
# .bim files
m2 <- BEDMatrix(path, n = 50, p = 1000)

# Alternatively, a BEDMatrix object can also be created using the `new`
# function
m3 <- new("BEDMatrix", path = path)
