#include "BEDMatrix_interface.h"

#include <R_ext/Rdynload.h>

/* compute_num_bytes_per_variant */

typedef int (*__compute_num_bytes_per_variant_funtype__)(
    int num_samples
);

int compute_num_bytes_per_variant(int num_samples) {
    static __compute_num_bytes_per_variant_funtype__ fun = NULL;
    if (fun == NULL) {
        fun = (__compute_num_bytes_per_variant_funtype__) R_GetCCallable("BEDMatrix", "compute_num_bytes_per_variant");
    }
    return fun(num_samples);
}


/* extract_genotype_linear */

typedef int (*__extract_genotype_linear_funtype__)(
    uint8_t *bed,
    ptrdiff_t k,
    int num_samples,
    int num_bytes_per_variant
);

int extract_genotype_linear(uint8_t *bed, ptrdiff_t k, int num_samples, int num_bytes_per_variant) {
    static __extract_genotype_linear_funtype__ fun = NULL;
    if (fun == NULL) {
        fun = (__extract_genotype_linear_funtype__) R_GetCCallable("BEDMatrix", "extract_genotype_linear");
    }
    return fun(bed, k, num_samples, num_bytes_per_variant);
}


/* extract_genotype_cartesian */

typedef int (*__extract_genotype_cartesian_funtype__)(
    uint8_t *bed,
    int i,
    int j,
    int num_samples
);

int extract_genotype_cartesian(uint8_t *bed, int i, int j, int num_bytes_per_variant) {
    static __extract_genotype_cartesian_funtype__ fun = NULL;
    if (fun == NULL) {
        fun = (__extract_genotype_cartesian_funtype__) R_GetCCallable("BEDMatrix", "extract_genotype_cartesian");
    }
    return fun(bed, i, j, num_bytes_per_variant);
}


/* recode_genotype */

typedef int (*__recode_genotype_funtype__)(
    int genotype,
    int na_value
);

int recode_genotype(int genotype, int na_value) {
    static __recode_genotype_funtype__ fun = NULL;
    if (fun == NULL) {
        fun = (__recode_genotype_funtype__) R_GetCCallable("BEDMatrix", "recode_genotype");
    }
    return fun(genotype, na_value);
}
