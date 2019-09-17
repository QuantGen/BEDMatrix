#ifndef BEDMATRIX_DEFINES_H
#define BEDMATRIX_DEFINES_H

#include <Rinternals.h>
#include <stdint.h>

struct BEDMatrix {
    int num_samples;
    int num_variants;
    uint8_t *data;
    size_t length;
};

struct mapped_region {
    void *addr;
    size_t length;
};

#endif
