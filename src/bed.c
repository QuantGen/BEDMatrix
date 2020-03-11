#include "bed.h"

#include "../inst/include/BEDMatrix.h"

int is_bed_file(uint8_t *bed) {
 // Check magic number
    if (!(bed[0] == 0x6c && bed[1] == 0x1b)) {
        return -1;
    }
 // Check mode: 00000001 indicates the default variant-major mode (i.e.
 // list all samples for first variant, all samples for second variant,
 // etc), 00000000 indicates the unsupported sample-major mode (i.e. list
 // all variants for the first sample, list all variants for the second
 // sample, etc)
    if (bed[2] != 0x01) {
        return -1;
    }
    return 0;
}

int has_valid_dimensions(size_t length, int num_samples, int num_variants) {
    int retval = 0;
 // File is a sequence of V blocks of N/4 (rounded up) bytes each, where V
 // is the number of variants and N is the number of samples.
    int num_bytes_per_variant = compute_num_bytes_per_variant(num_samples);
 // Check if given dimensions match the file
    if (((size_t) num_variants * num_bytes_per_variant) != (length - PLINK_BED_HEADER_LENGTH)) {
        retval = -1;
    }
    return retval;
}
