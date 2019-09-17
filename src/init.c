#include "BEDMatrix.h"

#include <R_ext/Rdynload.h>

static const R_CallMethodDef callEntries[] = {
    {"BEDMatrix_initialize", (DL_FUNC) &BEDMatrix_initialize, 3},
    {"BEDMatrix_extract_vector", (DL_FUNC) &BEDMatrix_extract_vector, 2},
    {"BEDMatrix_extract_matrix", (DL_FUNC) &BEDMatrix_extract_matrix, 3},
    {NULL, NULL, 0}
};

void R_init_BEDMatrix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, callEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
