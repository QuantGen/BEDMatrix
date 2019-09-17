#include "BEDMatrix_defines.h"

int compute_num_bytes_per_variant(
    int num_samples
);

int extract_genotype_linear(
    uint8_t *bed,
    ptrdiff_t k,
    int num_samples,
    int num_bytes_per_variant
);

int extract_genotype_cartesian(
    uint8_t *bed,
    int i,
    int j,
    int num_bytes_per_variant
);

int recode_genotype(
    int genotype,
    int na_value
);
