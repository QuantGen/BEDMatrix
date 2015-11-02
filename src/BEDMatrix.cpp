#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <string>

Rcpp::IntegerMatrix& preserve_dimnames(const Rcpp::List& x, Rcpp::IntegerMatrix& out, const Rcpp::IntegerVector& i, const Rcpp::IntegerVector& j) {
    Rcpp::List out_dimnames = Rcpp::List::create(
        R_NilValue,
        R_NilValue
    );
    Rcpp::List in_dimnames = x.attr("dnames");
    Rcpp::RObject in_rownames = in_dimnames[0];
    Rcpp::RObject in_colnames = in_dimnames[1];
    if (!in_rownames.isNULL()) {
        out_dimnames[0] = Rcpp::CharacterVector(in_rownames)[i];
    }
    if (!in_colnames.isNULL()) {
        out_dimnames[1] = Rcpp::CharacterVector(in_colnames)[j];
    }
    out.attr("dimnames") = out_dimnames;
    return out;
}

class BEDMatrix {
    public:
        BEDMatrix(std::string path, std::size_t n, std::size_t p);
        Rcpp::IntegerVector vector_subset(Rcpp::List x, Rcpp::IntegerVector i);
        Rcpp::IntegerMatrix matrix_subset(Rcpp::List x, Rcpp::IntegerVector i, Rcpp::IntegerVector j);
        std::size_t get_nrow();
        std::size_t get_ncol();
    private:
        BEDMatrix(const BEDMatrix&);
        BEDMatrix& operator=(const BEDMatrix&);
        int get_genotype(std::size_t i, std::size_t j);
        boost::iostreams::mapped_file_source file;
        const char* file_data;
        std::size_t nrow;
        std::size_t ncol;
        unsigned short int byte_padding; // Each new "row" starts a new byte.
        static const unsigned short int length_header;
};

BEDMatrix::BEDMatrix(std::string path, std::size_t n, std::size_t p) : file(path), nrow(n), ncol(p), byte_padding((n % 4 == 0) ? 0 : 4 - (n % 4)) {
    if (!this->file.is_open()) {
        Rcpp::stop("File not found.");
    }
    this->file_data = this->file.data();
    // Check magic number.
    if (!(this->file_data[0] == '\x6C' && this->file_data[1] == '\x1B')) {
        Rcpp::stop("File is not a binary PED file.");
    }
    // Check mode: 00000001 indicates the default SNP-major mode (i.e.
    // list all individuals for first SNP, all individuals for second
    // SNP, etc), 00000000 indicates the unsupported individual-major
    // mode (i.e. list all SNPs for the first individual, list all SNPs
    // for the second individual, etc).
    if (this->file_data[2] != '\x01') {
        Rcpp::stop("Individual-major mode is not supported.");
    }
    // Get number of bytes.
    const std::size_t num_bytes = this->file.size();
    // Check if given dimensions match the file.
    if ((this->nrow * this->ncol) + (this->byte_padding * this->ncol) != (num_bytes - this->length_header) * 4) {
        Rcpp::stop("n or p does not match the dimensions of the file.");
    }
}

int BEDMatrix::get_genotype(std::size_t i, std::size_t j) {
    // Reduce two-dimensional index to one-dimensional index with the mode.
    std::size_t which_pos = (j * this->nrow) + i + (this->byte_padding * j);
    // Every byte encodes 4 genotypes, find the one of interest.
    std::size_t which_byte = std::floor(which_pos / 4);
    // Find genotype in byte.
    unsigned short int which_genotype = (which_pos % 4) * 2;
    // Read in the whole byte.
    char genotypes = this->file_data[which_byte + this->length_header];
    // Remove the other genotypes by shifting the genotype of interest
    // to the end of the byte and masking with 00000011.
    char genotype = genotypes >> which_genotype & 3;
    // Remap genotype value.
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

Rcpp::IntegerVector BEDMatrix::vector_subset(Rcpp::List x, Rcpp::IntegerVector i) {
    // Check if index is out of bounds.
    if (Rcpp::is_true(Rcpp::any(i > this->nrow * this->ncol))) {
        Rcpp::stop("Invalid dimensions.");
    }
    // Convert from 1-index to 0-index.
    i = i - 1;
    // Keep size of i.
    std::size_t size_i = i.size();
    // Reserve output vector.
    Rcpp::IntegerVector out(size_i);
    // Iterate over indexes.
    for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {
        out(idx_i) = this->get_genotype(i[idx_i] % this->nrow, i[idx_i] / this->nrow);
    }
    return out;
}

Rcpp::IntegerMatrix BEDMatrix::matrix_subset(Rcpp::List x, Rcpp::IntegerVector i, Rcpp::IntegerVector j) {
    // Check if indexes are out of bounds.
    if (Rcpp::is_true(Rcpp::any(i > this->nrow)) || Rcpp::is_true(Rcpp::any(j > this->ncol))) {
        Rcpp::stop("Invalid dimensions.");
    }
    // Convert from 1-index to 0-index.
    i = i - 1;
    j = j - 1;
    // Keep sizes of i and j.
    std::size_t size_i = i.size();
    std::size_t size_j = j.size();
    // Reserve output matrix.
    Rcpp::IntegerMatrix out(size_i, size_j);
    preserve_dimnames(x, out, i, j);
    // Iterate over column indexes.
    for (std::size_t idx_j = 0; idx_j < size_j; idx_j++) {
        // Iterate over row indexes.
        for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {
            out(idx_i, idx_j) = this->get_genotype(i[idx_i], j[idx_j]);
        }
    }
    return out;
}


std::size_t BEDMatrix::get_nrow() {
    return this->nrow;
};

std::size_t BEDMatrix::get_ncol() {
    return this->ncol;
};

const unsigned short int BEDMatrix::length_header = 3;

RCPP_MODULE(mod_BEDMatrix) {

    using namespace Rcpp;

    class_<BEDMatrix>("BEDMatrix_")
    .constructor<std::string, std::size_t, std::size_t>()
    .method("vectorSubset", &BEDMatrix::vector_subset)
    .method("matrixSubset", &BEDMatrix::matrix_subset)
    .property("n", &BEDMatrix::get_nrow)
    .property("p", &BEDMatrix::get_ncol)
    ;

}
