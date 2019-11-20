#include <stddef.h>
#include <stdint.h>

int is_bed_file(
    uint8_t *bed
);

int has_valid_dimensions(
    size_t length,
    int num_samples,
    int num_variants
);
