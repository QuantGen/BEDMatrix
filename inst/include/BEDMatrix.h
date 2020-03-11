#ifndef BEDMATRIX_H
#define BEDMATRIX_H

#define PLINK_BED_HEADER_LENGTH 3
#define PLINK_BED_GENOTYPES_PER_BYTE 4

#include <math.h> // ceil
#include <stddef.h> // size_t, ptrdiff_t
#include <stdint.h> // uint8_t

#if defined(__GNUC__)
#define BEDMATRIX_EXPORT(declaration) __attribute__ ((unused)) static declaration
#else
#define BEDMATRIX_EXPORT(declaration) static declaration
#endif

struct BEDMatrix {
    int num_samples;
    int num_variants;
    uint8_t *data;
    size_t length;
};

BEDMATRIX_EXPORT(int compute_num_bytes_per_variant(int num_samples)) {
    return ceil((double) num_samples / PLINK_BED_GENOTYPES_PER_BYTE);
}

/**
 * Extract the genotype of the ith sample and the jth variant from a
 * variant-major .bed file
 * (https://www.cog-genomics.org/plink/1.9/formats#bed).
 *
 * The genotype is coded as follows: 00 indicates homozygous minor allele, 10
 * indicates heterozygous, 11 indicates homozygous major allele, and 01
 * indicates missing.
 *
 * 'num_bytes_per_variant' needs to be precomputed using the
 * 'compute_num_bytes_per_variant' function.
 */
BEDMATRIX_EXPORT(int extract_genotype_cartesian(uint8_t *bed,
                                                int i,
                                                int j,
                                                int num_bytes_per_variant)) {
 // Find corresponding byte
    int which_byte = i / PLINK_BED_GENOTYPES_PER_BYTE;
 // Find corresponding genotype position within byte
    int which_genotype = i % PLINK_BED_GENOTYPES_PER_BYTE;
 // Get byte from mapping
    uint8_t genotypes = bed[PLINK_BED_HEADER_LENGTH + ((ptrdiff_t) j * num_bytes_per_variant + which_byte)];
 // Extract genotype from byte by shifting bit pair of interest to the
 // right, then masking with 11
    return genotypes >> (2 * which_genotype) & 0x03;
}

BEDMATRIX_EXPORT(int extract_genotype_linear(uint8_t *bed,
                                             ptrdiff_t k,
                                             int num_samples,
                                             int num_bytes_per_variant)) {
    return extract_genotype_cartesian(
        bed,
        k % num_samples,
        k / num_samples,
        num_bytes_per_variant
    );
}

/**
 * Recode genotype value to count minor alleles, i.e., 0 indicates homozygous
 * major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
 * allele. A coding for the missing value needs to be provided in 'na_value'.
 * If used within R, 'NA_INTEGER' should be used.
 */
BEDMATRIX_EXPORT(int recode_genotype(int genotype, int na_value)) {
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
