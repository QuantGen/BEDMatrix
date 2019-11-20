#include <Rinternals.h>

SEXP BEDMatrix_initialize(
    SEXP path,
    SEXP n,
    SEXP p
);

SEXP BEDMatrix_extract_vector(
    SEXP xptr,
    SEXP k
);

SEXP BEDMatrix_extract_matrix(
    SEXP xptr,
    SEXP i,
    SEXP j
);
