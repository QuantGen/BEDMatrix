\docType{class}
\name{BEDMatrix-class}
\alias{BEDMatrix-class}
\alias{show,BEDMatrix-method}
\title{A Class to Extract Genotypes from a PLINK .bed File}
\description{
    \code{BEDMatrix} is a class that maps a
    \href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK .bed} file
    into memory and behaves similarly to a regular \code{matrix} by
    implementing key methods such as \code{[}, \code{dim}, and \code{dimnames}.
    Subsets are extracted directly and on-demand from the .bed file without
    loading the entire file into memory.
}
\details{
    The subsets extracted from a \code{BEDMatrix} object are coded similarly to
    \href{https://www.cog-genomics.org/plink2/formats#raw}{.raw} files
    (generated with the \code{--recodeA} argument in
    \href{https://www.cog-genomics.org/plink2/}{PLINK}): \code{0} indicates
    homozygous major allele, \code{1} indicates heterozygous, and \code{2}
    indicates homozygous minor allele.

    Internally, this class is an S4 class with the following slots that should
    not be relied upon in actual code: \code{xptr}, \code{dims}, \code{dnames},
    and \code{path}. The .bed file is mapped into memory using \code{mmap} on
    Unix and \code{MapViewOfFile} on Windows.
}
\section{Slots}{
    \describe{
        \item{\code{xptr}:}{
            An external pointer to the underlying \code{C} code.
        }
        \item{\code{dims}:}{
            An integer vector specifying the number of samples and variants as
            determined by the the accompanying
            \href{https://www.cog-genomics.org/plink2/formats#fam}{.fam} and
            \href{https://www.cog-genomics.org/plink2/formats#bim}{.bim} files
            or by the \code{n} and \code{p} parameters of the \code{BEDMatrix}
            constructor function.
        }
        \item{\code{dnames}:}{
            A list describing the row names and column names of the object as
            determined by the accompanying
            \href{https://www.cog-genomics.org/plink2/formats#fam}{.fam} and
            \href{https://www.cog-genomics.org/plink2/formats#bim}{.bim} files,
            or \code{NULL} if the \code{n} and \code{p} parameters of the
            \code{BEDMatrix} constructor function were provided.
        }
        \item{\code{path}:}{
            A character string containing the path to the .bed file.
        }
    }
}
\section{Methods}{
    \describe{
        \item{\code{[}:}{Extract parts of an object}
        \item{\code{dim}:}{Retrieve the dimension of an object}
        \item{\code{dimnames}:}{Retrieve the dimnames of an object}
        \item{\code{dimnames<-}:}{Set the dimnames of an object}
        \item{\code{as.matrix}:}{Turn the object into a matrix}
        \item{\code{is.matrix}:}{Test if the object is a matrix}
        \item{\code{length}:}{Get the length of an object}
        \item{\code{str}:}{Display the internal structure of an object}
        \item{\code{show}:}{Display the object}
    }
}
\seealso{
    \code{\link{BEDMatrix}} to create a \code{BEDMatrix} object from a .bed
    file, \code{\link{BEDMatrix-package}} to learn more about the
    \code{BEDMatrix} package,
    \code{\link[LinkedMatrix:LinkedMatrix-package]{LinkedMatrix}} to link
    several \code{BEDMatrix} objects together.
}
\examples{
# Get the path to the example .bed file
path <- system.file("extdata", "example.bed",
                    package = "BEDMatrix")

# Create a BEDMatrix object the example .bed file
m <- BEDMatrix(path)

# Get the dimensions of the BEDMatrix object
dim(m)

# Get the row names of the BEDMatrix object
rownames(m)

# Get the column names of the BEDMatrix object
colnames(m)

# Extract genotypes for the specified sample(s)
m[1, ]
m[1:3, ]
m["per0_per0", ]
m[c("per0_per0", "per1_per1", "per2_per2"), ]

# Extract genotypes for a particular variant
m[, 1]
m[, c("snp0_A", "snp1_C", "snp2_G")]

# Extract genotypes for the specified samples and variants
m[
    c("per0_per0", "per1_per1", "per2_per2"),
    c("snp0_A", "snp1_C", "snp2_G")
]
}
\keyword{classes}
