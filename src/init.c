#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP BEDMatrix_new(SEXP, SEXP, SEXP);
extern SEXP BEDMatrix_extract_matrix(SEXP, SEXP, SEXP);
extern SEXP BEDMatrix_extract_vector(SEXP, SEXP);

static const R_CallMethodDef callEntries[] = {
    {"BEDMatrix_new", (DL_FUNC) &BEDMatrix_new, 3},
    {"BEDMatrix_extract_matrix", (DL_FUNC) &BEDMatrix_extract_matrix, 3},
    {"BEDMatrix_extract_vector", (DL_FUNC) &BEDMatrix_extract_vector, 2},
    {NULL, NULL, 0}
};

void R_init_BEDMatrix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
