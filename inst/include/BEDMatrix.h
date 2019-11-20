#ifndef BEDMATRIX_H
#define BEDMATRIX_H

#define PLINK_BED_HEADER_LENGTH 3
#define PLINK_BED_GENOTYPES_PER_BYTE 4

#include <math.h> // for ceil
#include <stddef.h>
#include <stdint.h>

struct BEDMatrix {
    int num_samples;
    int num_variants;
    uint8_t *data;
    size_t length;
};

static int compute_num_bytes_per_variant(int num_samples) {
    return ceil((double) num_samples / PLINK_BED_GENOTYPES_PER_BYTE);
};

static int extract_genotype_cartesian(uint8_t *bed, int i, int j, int num_bytes_per_variant) {
 // Each byte contains 4 genotypes; adjust indices
    int which_byte = i / PLINK_BED_GENOTYPES_PER_BYTE;
    int which_genotype = i % PLINK_BED_GENOTYPES_PER_BYTE;
 // Get byte from mapping
    uint8_t genotypes = bed[PLINK_BED_HEADER_LENGTH + ((ptrdiff_t) j * num_bytes_per_variant + which_byte)];
 // Extract genotype from byte by shifting bit pair of interest to the
 // right, then masking with 11
    return genotypes >> (2 * which_genotype) & 0x03;
}

static int extract_genotype_linear(uint8_t *bed, ptrdiff_t k, int num_samples, int num_bytes_per_variant) {
    return extract_genotype_cartesian(
        bed,
        k % num_samples,
        k / num_samples,
        num_bytes_per_variant
    );
}

static int recode_genotype(int genotype, int na_value) {
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

#endif
