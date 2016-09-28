# Create an example BEDMatrix object
m <- BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))

# Get the dimensions of the example BEDMatrix object
dim(m)

# Extract a subset of the example BEDMatrix object
m[1:3, ]
