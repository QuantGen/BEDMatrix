#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP C_extract_matrix(SEXP, SEXP, SEXP);
extern SEXP C_extract_vector(SEXP, SEXP);
extern SEXP C_new(SEXP, SEXP, SEXP);

static const R_CallMethodDef callEntries[] = {
    {"C_extract_matrix", (DL_FUNC) &C_extract_matrix, 3},
    {"C_extract_vector", (DL_FUNC) &C_extract_vector, 2},
    {"C_new", (DL_FUNC) &C_new, 3},
    {NULL, NULL, 0}
};

void R_init_BEDMatrix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
