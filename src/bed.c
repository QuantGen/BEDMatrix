#include "BEDMatrix.h"

#include <errno.h>
#include <math.h>

#define PLINK_BED_HEADER_LENGTH 3
#define PLINK_BED_GENOTYPES_PER_BYTE 4

int compute_num_bytes_per_variant(int num_samples) {
    return ceil((double) num_samples / PLINK_BED_GENOTYPES_PER_BYTE);
}

int is_bed_file(uint8_t *bed) {
 // Check magic number
    if (!(bed[0] == 0x6c && bed[1] == 0x1b)) {
        errno = 1;
        return -1;
    }
 // Check mode: 00000001 indicates the default variant-major mode (i.e.
 // list all samples for first variant, all samples for second variant,
 // etc), 00000000 indicates the unsupported sample-major mode (i.e. list
 // all variants for the first sample, list all variants for the second
 // sample, etc)
    if (bed[2] != 0x01) {
        errno = 2;
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

int extract_genotype_linear(uint8_t *bed, ptrdiff_t k, int num_samples, int num_bytes_per_variant) {
    return extract_genotype_cartesian(
        bed,
        k % num_samples,
        k / num_samples,
        num_bytes_per_variant
    );
}

int extract_genotype_cartesian(uint8_t *bed, int i, int j, int num_bytes_per_variant) {
 // Each byte contains 4 genotypes; adjust indices
    int which_byte = i / PLINK_BED_GENOTYPES_PER_BYTE;
    int which_genotype = i % PLINK_BED_GENOTYPES_PER_BYTE;
 // Get byte from mapping
    uint8_t genotypes = bed[PLINK_BED_HEADER_LENGTH + ((ptrdiff_t) j * num_bytes_per_variant + which_byte)];
 // Extract genotype from byte by shifting bit pair of interest to the
 // right, then masking with 11
    return genotypes >> (2 * which_genotype) & 0x03;
}

int recode_genotype(int genotype, int na_value) {
 // Recode genotype value to resemble RAW file, i.e. 0 indicates homozygous
 // major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
 // allele. In BED, the coding is different: homozygous minor allele is 0 (00)
 // and homozygous major allele is 3 (11). Each byte is read backwards.
    int coding = na_value; // missing
    if (genotype == 0) {
        coding = 2; // homozygous AA
    } else if (genotype == 3) {
        coding = 0; // homozygous BB
    } else if (genotype == 2) {
        coding = 1; // heterozygous AB
    }
    return coding;
}
