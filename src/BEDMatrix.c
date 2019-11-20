#include "BEDMatrix.h"

#include "bed.h"
#include "mapping.h"

#include "../inst/include/BEDMatrix.h"

#include <R_ext/RS.h>
#include <R_ext/Utils.h>

#define INTERRUPT_INTERVAL 10000000

static void BEDMatrix_finalize(SEXP xptr) {
    struct BEDMatrix *state = R_ExternalPtrAddr(xptr);
    if (state == NULL) {
        return;
    }
    struct mapped_region mapped_file;
    mapped_file.addr = state->data;
    mapped_file.length = state->length;
    unmap_file(&mapped_file); // ignore errors
    Free(state);
    R_ClearExternalPtr(xptr);
}

SEXP BEDMatrix_initialize(SEXP path, SEXP n, SEXP p) {
    const char *expanded_filename = R_ExpandFileName(CHAR(Rf_asChar(path)));
    int nrows = Rf_asInteger(n);
    int ncols = Rf_asInteger(p);
 // Map file
    struct mapped_region mapped_file;
    if (map_file(expanded_filename, &mapped_file, 'r') == -1) {
        Rf_error("Could not map file.");
    }
 // Test if file is a valid .bed file
    if (is_bed_file(mapped_file.addr) == -1) {
        unmap_file(&mapped_file); // ignore errors
        Rf_error("File is not a PLINK .bed file.");
    }
 // Test if n and p correspond to length
    if (has_valid_dimensions(mapped_file.length, nrows, ncols) == -1) {
        unmap_file(&mapped_file); // ignore errors
        Rf_error("n or p does not match the dimensions of the file.");
    }
 // Create state
    struct BEDMatrix *state = Calloc(1, struct BEDMatrix);
    state->num_samples = nrows;
    state->num_variants = ncols;
    state->data = mapped_file.addr;
    state->length = mapped_file.length;
 // Create external pointer
    SEXP xptr = PROTECT(R_MakeExternalPtr(
        state,
        R_NilValue,
        R_NilValue
    ));
 // Register finalizer
    R_RegisterCFinalizerEx(xptr, BEDMatrix_finalize, 1);
    UNPROTECT(1); // xptr
    return xptr;
}

SEXP BEDMatrix_extract_vector(SEXP xptr, SEXP k) {
    struct BEDMatrix *state = R_ExternalPtrAddr(xptr);
    if (state == NULL) {
        Rf_error("BEDMatrix instance has been unmapped.");
    }
    ptrdiff_t nx = (ptrdiff_t) state->num_samples * state->num_variants;
    int num_bytes_per_variant = compute_num_bytes_per_variant(state->num_samples);
    R_xlen_t nk = Rf_xlength(k);
    SEXP out = PROTECT(Rf_allocVector(INTSXP, nk));
    int *pout = INTEGER(out);
    // index can be either integer or double, depending on whether an index
    // value is greater than the largest integer
    if (TYPEOF(k) == INTSXP) {
        int *pk = INTEGER(k);
        for (ptrdiff_t ck = 0; ck < nk; ck++) {
            int kk = pk[ck];
            if (0 < kk && kk <= nx) {
                kk--;
                pout[ck] = recode_genotype(extract_genotype_linear(
                    state->data,
                    kk,
                    state->num_samples,
                    num_bytes_per_variant
                ), NA_INTEGER);
            } else {
                pout[ck] = NA_INTEGER;
            }
            if (ck % INTERRUPT_INTERVAL == 0) {
                R_CheckUserInterrupt();
            }
        }
    } else {
        double *pk = REAL(k);
        for (ptrdiff_t ck = 0; ck < nk; ck++) {
            double dk = pk[ck];
            ptrdiff_t kk = (ptrdiff_t) (dk - 1);
            if (R_FINITE(dk) && 0 <= kk && kk < nx) {
                pout[ck] = recode_genotype(extract_genotype_linear(
                    state->data,
                    kk,
                    state->num_samples,
                    num_bytes_per_variant
                ), NA_INTEGER);
            } else {
                pout[ck] = NA_INTEGER;
            }
            if (ck % INTERRUPT_INTERVAL == 0) {
                R_CheckUserInterrupt();
            }
        }
    }
    UNPROTECT(1); // out
    return out;
}

SEXP BEDMatrix_extract_matrix(SEXP xptr, SEXP i, SEXP j) {
    struct BEDMatrix *state = R_ExternalPtrAddr(xptr);
    if (state == NULL) {
        Rf_error("BEDMatrix instance has been unmapped.");
    }
    int num_bytes_per_variant = compute_num_bytes_per_variant(state->num_samples);
    R_xlen_t ni = Rf_xlength(i);
    int *pi = INTEGER(i);
    R_xlen_t nj = Rf_xlength(j);
    int *pj = INTEGER(j);
    SEXP out = PROTECT(Rf_allocMatrix(INTSXP, ni, nj));
    int *pout = INTEGER(out);
    ptrdiff_t cint = 0;
    for (int cj = 0; cj < nj; cj++) {
        int jj = pj[cj];
        for (int ci = 0; ci < ni; ci++) {
            int ii = pi[ci];
            if (ii == NA_INTEGER || jj == NA_INTEGER) {
                pout[(ptrdiff_t) cj * ni + ci] = NA_INTEGER;
            } else {
                int genotype = recode_genotype(extract_genotype_cartesian(
                    state->data,
                    ii - 1,
                    jj - 1,
                    num_bytes_per_variant
                ), NA_INTEGER);
                pout[(ptrdiff_t) cj * ni + ci] = genotype;
            }
            if (cint % INTERRUPT_INTERVAL == 0) {
                R_CheckUserInterrupt();
            }
            cint++;
        }
    }
    UNPROTECT(1); // out
    return out;
}
