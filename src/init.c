#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP BEDMatrix__extract_matrix(SEXP, SEXP, SEXP);
extern SEXP BEDMatrix__extract_vector(SEXP, SEXP);
extern SEXP BEDMatrix__new(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"BEDMatrix__extract_matrix", (DL_FUNC) &BEDMatrix__extract_matrix, 3},
    {"BEDMatrix__extract_vector", (DL_FUNC) &BEDMatrix__extract_vector, 2},
    {"BEDMatrix__new",            (DL_FUNC) &BEDMatrix__new,            3},
    {NULL, NULL, 0}
};

void R_init_BEDMatrix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
