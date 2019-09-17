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
    R_RegisterCCallable("BEDMatrix", "compute_num_bytes_per_variant", (DL_FUNC) &compute_num_bytes_per_variant);
    R_RegisterCCallable("BEDMatrix", "extract_genotype_linear", (DL_FUNC) &extract_genotype_linear);
    R_RegisterCCallable("BEDMatrix", "extract_genotype_cartesian", (DL_FUNC) &extract_genotype_cartesian);
    R_RegisterCCallable("BEDMatrix", "recode_genotype", (DL_FUNC) &recode_genotype);
}
