#ifndef BEDMATRIX_H
#define BEDMATRIX_H

#include <Rinternals.h>
#include <stdint.h>

/* bed.c */

int compute_num_bytes_per_variant(
    int num_samples
);

int is_bed_file(
    uint8_t *bed
);

int has_valid_dimensions(
    size_t length,
    int num_samples,
    int num_variants
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


/* mapping.c */

struct mapped_region {
    void *addr;
    size_t length;
};

int map_file(const char *pathname, struct mapped_region *mapped_region, char mode);

int unmap_file(struct mapped_region *mapped_region);


/* BEDMatrix.c */

struct BEDMatrix {
    int num_samples;
    int num_variants;
    uint8_t *data;
    size_t length;
};

SEXP BEDMatrix_new(
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

#endif
