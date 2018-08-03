// [[Rcpp::depends(BH)]]

#define PLINK_BED_HEADER_LENGTH 3
#define PLINK_BED_GENOTYPES_PER_BYTE 4

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/exceptions.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <string>

class BEDMatrix {
    public:
        BEDMatrix(std::string path, std::size_t n, std::size_t p);
        Rcpp::IntegerVector extract_vector(Rcpp::IntegerVector i);
        Rcpp::IntegerMatrix extract_matrix(Rcpp::IntegerVector i, Rcpp::IntegerVector j);
    private:
        BEDMatrix(const BEDMatrix&);
        BEDMatrix& operator=(const BEDMatrix&);
        int get_genotype(std::size_t i, std::size_t j);
        boost::interprocess::file_mapping file;
        boost::interprocess::mapped_region file_region;
        uint8_t* file_data;
        std::size_t nrow;
        std::size_t ncol;
};

BEDMatrix::BEDMatrix(std::string path, std::size_t n, std::size_t p) : nrow(n), ncol(p) {
    try {
        this->file = boost::interprocess::file_mapping(path.c_str(), boost::interprocess::read_only);
    } catch(const boost::interprocess::interprocess_exception& e) {
        throw std::runtime_error("File not found.");
    }
    this->file_region = boost::interprocess::mapped_region(this->file, boost::interprocess::read_only);
    this->file_data = static_cast<uint8_t*>(this->file_region.get_address());
    // Check magic number
    if (!(this->file_data[0] == '\x6C' && this->file_data[1] == '\x1B')) {
        throw std::runtime_error("File is not a binary PED file.");
    }
    // Check mode: 00000001 indicates the default variant-major mode (i.e.
    // list all samples for first variant, all samples for second variant,
    // etc), 00000000 indicates the unsupported sample-major mode (i.e. list
    // all variants for the first sample, list all variants for the second
    // sample, etc)
    if (this->file_data[2] != '\x01') {
        throw std::runtime_error("Sample-major mode is not supported.");
    }
    // Get number of bytes
    const std::size_t num_bytes = this->file_region.get_size();
    // Check if given dimensions match the file
    if ((this->ncol * ceil((double) this->nrow / PLINK_BED_GENOTYPES_PER_BYTE)) != (num_bytes - PLINK_BED_HEADER_LENGTH)) {
        throw std::runtime_error("n or p does not match the dimensions of the file.");
    }
}

int BEDMatrix::get_genotype(std::size_t i, std::size_t j) {
    // Each byte encodes 4 genotypes; adjust indices
    std::size_t i_bytes = i / PLINK_BED_GENOTYPES_PER_BYTE;
    std::size_t n_bytes = this->nrow / PLINK_BED_GENOTYPES_PER_BYTE + (this->nrow % PLINK_BED_GENOTYPES_PER_BYTE != 0); // fast ceil for int
    std::size_t i_genotypes = 2 * (i - i_bytes * PLINK_BED_GENOTYPES_PER_BYTE);
    // Load byte from map
    uint8_t genotypes = this->file_data[PLINK_BED_HEADER_LENGTH + (j * n_bytes + i_bytes)];
    // Extract genotypes from byte by shifting the genotype of interest to the
    // end of the byte and masking with 00000011
    uint8_t genotype = genotypes >> i_genotypes & 3;
    // Remap genotype value to resemble RAW file, i.e. 0 indicates homozygous
    // major allele, 1 indicates heterozygous, and 2 indicates homozygous minor
    // allele. In BED, the coding is different: homozygous minor allele is 0
    // (00) and homozygous major allele is 3 (11). Each byte is read backwards.
    int mapping = NA_INTEGER; // missing
    if (genotype == 0) {
        mapping = 2; // homozygous AA
    } else if (genotype == 3) {
        mapping = 0; // homozygous BB
    } else if (genotype == 2) {
        mapping = 1; // heterozygous AB
    }
    return mapping;
}

Rcpp::IntegerVector BEDMatrix::extract_vector(Rcpp::IntegerVector i) {
    // Convert from 1-index to 0-index
    Rcpp::IntegerVector i0(i - 1);
    // Keep size of i
    std::size_t size_i = i.size();
    // Reserve output vector
    Rcpp::IntegerVector out(size_i);
    // Get bounds
    std::size_t bounds = this->nrow * this->ncol;
    // Iterate over indexes
    for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {
        if (Rcpp::IntegerVector::is_na(i0[idx_i]) || static_cast<std::size_t>(i0[idx_i]) >= bounds) {
            out(idx_i) = NA_INTEGER;
        } else {
            out(idx_i) = this->get_genotype(i0[idx_i] % this->nrow, i0[idx_i] / this->nrow);
        }
    }
    return out;
}

Rcpp::IntegerMatrix BEDMatrix::extract_matrix(Rcpp::IntegerVector i, Rcpp::IntegerVector j) {
    // Check if indexes are out of bounds
    if (Rcpp::is_true(Rcpp::any(i > this->nrow)) || Rcpp::is_true(Rcpp::any(j > this->ncol))) {
        throw std::runtime_error("subscript out of bounds");
    }
    // Convert from 1-index to 0-index
    Rcpp::IntegerVector i0(i - 1);
    Rcpp::IntegerVector j0(j - 1);
    // Keep sizes of i and j
    std::size_t size_i = i.size();
    std::size_t size_j = j.size();
    // Reserve output matrix
    Rcpp::IntegerMatrix out(size_i, size_j);
    // Iterate over column indexes
    for (std::size_t idx_j = 0; idx_j < size_j; idx_j++) {
        // Iterate over row indexes
        for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {
            if (Rcpp::IntegerVector::is_na(i0[idx_i]) || Rcpp::IntegerVector::is_na(j0[idx_j])) {
                out(idx_i, idx_j) = NA_INTEGER;
            } else {
                out(idx_i, idx_j) = this->get_genotype(i0[idx_i], j0[idx_j]);
            }
        }
    }
    return out;
}

// Export BEDMatrix::BEDMatrix
RcppExport SEXP C_new(SEXP path_, SEXP n_, SEXP p_) {
    // Convert inputs to appropriate C++ types
    std::string path = Rcpp::as<std::string>(path_);
    std::size_t n = Rcpp::as<std::size_t>(n_);
    std::size_t p = Rcpp::as<std::size_t>(p_);
    try {
        // Create a pointer to a BEDMatrix object and wrap it as an external
        // pointer
        Rcpp::XPtr<BEDMatrix> ptr(new BEDMatrix(path, n, p), true);
        // Return the external pointer to the R side
        return ptr;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
};

// Export BEDMatrix::extract_vector
RcppExport SEXP C_extract_vector(SEXP xp_, SEXP i_) {
    // Convert inputs to appropriate C++ types
    Rcpp::XPtr<BEDMatrix> ptr(xp_);
    Rcpp::IntegerVector i = Rcpp::as<Rcpp::IntegerVector>(i_);
    try {
        // Invoke the extract_vector function
        Rcpp::IntegerVector res = ptr->extract_vector(i);
        return res;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
};

// Export BEDMatrix::extract_matrix
RcppExport SEXP C_extract_matrix(SEXP xp_, SEXP i_, SEXP j_) {
    // Convert inputs to appropriate C++ types
    Rcpp::XPtr<BEDMatrix> ptr(xp_);
    Rcpp::IntegerVector i = Rcpp::as<Rcpp::IntegerVector>(i_);
    Rcpp::IntegerVector j = Rcpp::as<Rcpp::IntegerVector>(j_);
    try {
        // Invoke the extract_matrix function
        Rcpp::IntegerMatrix res = ptr->extract_matrix(i, j);
        return res;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
        return 0;
    }
};
