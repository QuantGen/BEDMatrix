#include <fstream>
#include <iostream>
#include <string>
#include <Rcpp.h>

Rcpp::IntegerMatrix& preserveDimnames(const Rcpp::List& x, Rcpp::IntegerMatrix& out, const Rcpp::IntegerVector& i, const Rcpp::IntegerVector& j) {
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
        BEDMatrix(std::string path, unsigned int n, unsigned int p);
        ~BEDMatrix();
        Rcpp::IntegerVector vector_subset(Rcpp::List x, Rcpp::IntegerVector i);
        Rcpp::IntegerMatrix matrix_subset(Rcpp::List x, Rcpp::IntegerVector i, Rcpp::IntegerVector j);
        unsigned int get_nrow();
        unsigned int get_ncol();
    private:
        BEDMatrix(const BEDMatrix&);
        BEDMatrix& operator=(const BEDMatrix&);
        int get_genotype(unsigned int i, unsigned int j);
        std::ifstream infile;
        unsigned int nrow;
        unsigned int ncol;
        unsigned int byte_padding; // Each new "row" starts a new byte.
        static const unsigned int length_header;
};

BEDMatrix::BEDMatrix(std::string path, unsigned int n, unsigned int p) : infile(path.c_str(), std::ios::binary), nrow(n), ncol(p), byte_padding((n % 4 == 0) ? 0 : 4 - (n % 4)) {
    if (!this->infile) {
        Rcpp::stop("File not found.");
    }
    char header[this->length_header];
    infile.read(header, this->length_header);
    // Check magic number.
    if (!(header[0] == '\x6C' && header[1] == '\x1B')) {
        Rcpp::stop("File is not a binary PED file.");
    }
    // Check mode: 00000001 indicates the default SNP-major mode (i.e.
    // list all individuals for first SNP, all individuals for second
    // SNP, etc), 00000000 indicates the unsupported individual-major
    // mode (i.e. list all SNPs for the first individual, list all SNPs
    // for the second individual, etc).
    if (header[2] != '\x01') {
        Rcpp::stop("Individual-major mode is not supported.");
    }
    // Get number of bytes.
    this->infile.seekg(0, infile.end);
    unsigned int num_bytes = infile.tellg();
    // Check if given dimensions match the file.
    if ((this->nrow * this->ncol) + (this->byte_padding * this->ncol) != (num_bytes - this->length_header) * 4) {
        Rcpp::stop("n or p does not match the dimensions of the file.");
    }
}

BEDMatrix::~BEDMatrix() {
    this->infile.close();
}

int BEDMatrix::get_genotype(unsigned int i, unsigned int j) {
    // Reduce two-dimensional index to one-dimensional index with the mode.
    unsigned int which_pos = (j * this->nrow) + i + (this->byte_padding * j);
    // Every byte encodes 4 genotypes, find the one of interest.
    unsigned int which_byte = std::floor(which_pos / 4);
    // Find genotype in byte.
    unsigned int which_genotype = (which_pos % 4) * 2;
    // Read in the whole byte.
    infile.seekg(which_byte + this->length_header);
    char genotypes = 0;
    infile.read(&genotypes, 1);
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
    unsigned int size_i = i.size();
    // Reserve output vector.
    Rcpp::IntegerVector out(size_i);
    // Iterate over indexes.
    for (unsigned int idx_i = 0; idx_i < size_i; idx_i++) {
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
    unsigned int size_i = i.size();
    unsigned int size_j = j.size();
    // Reserve output matrix.
    Rcpp::IntegerMatrix out(size_i, size_j);
    preserveDimnames(x, out, i, j);
    // Iterate over column indexes.
    for (unsigned int idx_j = 0; idx_j < size_j; idx_j++) {
        // Iterate over row indexes.
        for (unsigned int idx_i = 0; idx_i < size_i; idx_i++) {
            out(idx_i, idx_j) = this->get_genotype(i[idx_i], j[idx_j]);
        }
    }
    return out;
}


unsigned int BEDMatrix::get_nrow() {
    return this->nrow;
};

unsigned int BEDMatrix::get_ncol() {
    return this->ncol;
};

const unsigned int BEDMatrix::length_header = 3;

RCPP_MODULE(mod_BEDMatrix) {

    using namespace Rcpp;

    class_<BEDMatrix>("BEDMatrix_")
    .constructor<std::string, unsigned int, unsigned int>()
    .method("vectorSubset", &BEDMatrix::vector_subset)
    .method("matrixSubset", &BEDMatrix::matrix_subset)
    .property("n", &BEDMatrix::get_nrow)
    .property("p", &BEDMatrix::get_ncol)
    ;

}
