// [[Rcpp::depends(BH)]]

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/exceptions.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <string>

static const int plink_bed_header_length = 3;
static const int plink_bed_genotypes_per_byte = 4;

class BEDMatrix {
    public:
        BEDMatrix(std::string path, std::size_t n, std::size_t p);
        Rcpp::IntegerVector extract_vector(Rcpp::IntegerVector i);
        Rcpp::IntegerVector extract_vector(Rcpp::NumericVector i);
        Rcpp::IntegerMatrix extract_matrix(Rcpp::IntegerVector i, Rcpp::IntegerVector j);
    private:
        BEDMatrix(const BEDMatrix&);
        BEDMatrix& operator=(const BEDMatrix&);
        int get_genotype(std::size_t i, std::size_t j);
        boost::interprocess::file_mapping file;
        boost::interprocess::mapped_region file_region;
        uint8_t* file_data;
        std::size_t num_samples;
        std::size_t num_variants;
        std::size_t num_bytes_per_variant;
};

BEDMatrix::BEDMatrix(std::string path, std::size_t n, std::size_t p) : num_samples(n), num_variants(p) {
    try {
        this->file = boost::interprocess::file_mapping(path.c_str(), boost::interprocess::read_only);
    } catch(const boost::interprocess::interprocess_exception& e) {
        throw std::runtime_error("File not found.");
    }
    this->file_region = boost::interprocess::mapped_region(this->file, boost::interprocess::read_only);
    this->file_data = static_cast<uint8_t*>(this->file_region.get_address());
    // Check magic number
    if (!(this->file_data[0] == 0x6C && this->file_data[1] == 0x1B)) {
        throw std::runtime_error("File is not a PLINK .bed file.");
    }
    // Check mode: 00000001 indicates the default variant-major mode (i.e.
    // list all samples for first variant, all samples for second variant,
    // etc), 00000000 indicates the unsupported sample-major mode (i.e. list
    // all variants for the first sample, list all variants for the second
    // sample, etc)
    if (this->file_data[2] != 0x01) {
        throw std::runtime_error("Sample-major mode is not supported.");
    }
    // Get number of bytes
    const std::size_t num_bytes = this->file_region.get_size();
    // Determine the number of bytes per variant
    this->num_bytes_per_variant = ceil((this->num_samples + 0.0) / plink_bed_genotypes_per_byte);
    // Check if given dimensions match the file
    if ((this->num_variants * this->num_bytes_per_variant) != (num_bytes - plink_bed_header_length)) {
        throw std::runtime_error("n or p does not match the dimensions of the file.");
    }
}

int BEDMatrix::get_genotype(std::size_t i, std::size_t j) {
    // Each byte encodes 4 genotypes; determine which byte contains the
    // genotype at index [i, j]
    std::size_t which_byte = i / plink_bed_genotypes_per_byte;
    // Extract byte from memory-mapped region
    uint8_t genotypes = this->file_data[plink_bed_header_length + (j * this->num_bytes_per_variant + which_byte)];
    // Determine which genotype in the byte to extract
    std::size_t which_genotype = i % plink_bed_genotypes_per_byte;
    // Extract genotype from byte by shifting the genotype of interest to the
    // end of the byte and masking with 00000011
    uint8_t genotype = genotypes >> (2 * which_genotype) & 0x03;
    // Remap genotype value to resemble RAW file, i.e. 0 indicates homozygous
    // major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
    // allele. In BED, the coding is different: homozygous minor allele is 0
    // (00) and homozygous major allele is 3 (11). Each byte is read backwards.
    int mapping = NA_INTEGER; // missing
    if (genotype == 0x00) {
        mapping = 2; // homozygous AA
    } else if (genotype == 0x03) {
        mapping = 0; // homozygous BB
    } else if (genotype == 0x02) {
        mapping = 1; // heterozygous AB
    }
    return mapping;
}

Rcpp::IntegerVector BEDMatrix::extract_vector(Rcpp::IntegerVector i) {
    // Keep size of i
    R_xlen_t size_i = i.size();
    // Reserve output vector
    Rcpp::IntegerVector out(size_i);
    // Compute length
    R_xlen_t length = this->num_samples * this->num_variants;
    // Iterate over index
    for (R_xlen_t idx_i = 0; idx_i < size_i; idx_i++) {
        R_xlen_t ii = i[idx_i];
        if (0 < ii && ii <= length) {
            ii--;
            out(idx_i) = this->get_genotype(ii % this->num_samples, ii / this->num_samples);
        } else {
            out(idx_i) = NA_INTEGER;
        }
    }
    return out;
}

Rcpp::IntegerVector BEDMatrix::extract_vector(Rcpp::NumericVector i) {
    // Keep size of i
    R_xlen_t size_i = i.size();
    // Reserve output vector
    Rcpp::IntegerVector out(size_i);
    // Compute length
    R_xlen_t length = this->num_samples * this->num_variants;
    // Iterate over index
    for (R_xlen_t idx_i = 0; idx_i < size_i; idx_i++) {
        double di = i[idx_i];
        R_xlen_t ii = static_cast<R_xlen_t>(di - 1);
        if (R_FINITE(di) && 0 <= ii && ii < length) {
            out(idx_i) = this->get_genotype(ii % this->num_samples, ii / this->num_samples);
        } else {
            out(idx_i) = NA_INTEGER;
        }
    }
    return out;
}

/**
 * extract_matrix expects that i and j have been bound checked.
 */
Rcpp::IntegerMatrix BEDMatrix::extract_matrix(Rcpp::IntegerVector i, Rcpp::IntegerVector j) {
    // Keep sizes of i and j
    R_xlen_t size_i = i.size();
    R_xlen_t size_j = j.size();
    // Reserve output matrix
    Rcpp::IntegerMatrix out(size_i, size_j);
    // Iterate over column index
    for (R_xlen_t idx_j = 0; idx_j < size_j; idx_j++) {
        R_xlen_t jj = j[idx_j];
        // Iterate over row index
        for (R_xlen_t idx_i = 0; idx_i < size_i; idx_i++) {
            R_xlen_t ii = i[idx_i];
            if (ii == NA_INTEGER || jj == NA_INTEGER) {
                out(idx_i, idx_j) = NA_INTEGER;
            } else {
                out(idx_i, idx_j) = this->get_genotype(ii - 1, jj - 1);
            }
        }
    }
    return out;
}

/**
 * Export BEDMatrix::BEDMatrix
 */
RcppExport SEXP C_new(SEXP path_, SEXP n_, SEXP p_) {
    std::string path = Rcpp::as<std::string>(path_);
    std::size_t n = Rcpp::as<std::size_t>(n_);
    std::size_t p = Rcpp::as<std::size_t>(p_);
    try {
        // Create a BEDMatrix object and return its pointer to the R side as an
        // external pointer
        Rcpp::XPtr<BEDMatrix> ptr(new BEDMatrix(path, n, p), true);
        return ptr;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
}

/**
 * Export BEDMatrix::extract_vector
 */
RcppExport SEXP C_extract_vector(SEXP xp_, SEXP i_) {
    Rcpp::XPtr<BEDMatrix> ptr(xp_);
    try {
        Rcpp::IntegerVector res;
        // index can be either integer or double, depending on whether an index
        // value is greater than the largest integer
        if (TYPEOF(i_) == INTSXP) {
            res = ptr->extract_vector(Rcpp::IntegerVector(i_));
        } else {
            res = ptr->extract_vector(Rcpp::NumericVector(i_));
        }
        return res;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
}

/**
 * Export BEDMatrix::extract_matrix
 */
RcppExport SEXP C_extract_matrix(SEXP xp_, SEXP i_, SEXP j_) {
    Rcpp::XPtr<BEDMatrix> ptr(xp_);
    Rcpp::IntegerVector i(i_);
    Rcpp::IntegerVector j(j_);
    try {
        Rcpp::IntegerMatrix res = ptr->extract_matrix(i, j);
        return res;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
}
