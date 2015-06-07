#include <Rcpp.h>

#include <string>
#include <iostream>
#include <fstream>

// [[Rcpp::export]]
Rcpp::IntegerMatrix subsetBED(Rcpp::List x, Rcpp::IntegerVector i, Rcpp::IntegerVector j) {
  // Check if x is a BEDMatrix.
  if (!x.inherits("BEDMatrix")) {
    Rcpp::stop("x must be a BEDMatrix.");
  }
  std::string path = x.attr("path");
  int n = x.attr("n");
  int p = x.attr("p");
  // Check if indexes are out of bounds.
  if (Rcpp::is_true(Rcpp::any(i > n)) || Rcpp::is_true(Rcpp::any(j > p))) {
    Rcpp::stop("Invalid dimensions.");
  }
  // Convert from 1-index to 0-index.
  i = i - 1;
  j = j - 1;
  // Keep sizes of i and j.
  int size_i = i.size();
  int size_j = j.size();
  // Reserve output matrix.
  Rcpp::IntegerMatrix out (size_i, size_j);
  // Preserve dimnames.
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
  // Open BED file.
  std::ifstream infile (path.c_str(), std::ios::binary);
  if (infile) {
    char *header = new char[3];
    infile.read(header, 3);
    // Check magic number.
    if (header[0] == '\x6C' && header[1] == '\x1B') {
      // Check mode: 00000001 indicates the default SNP-major mode (i.e.
      // list all individuals for first SNP, all individuals for second
      // SNP, etc), 00000000 indicates the unsupported individual-major
      // mode (i.e. list all SNPs for the first individual, list all SNPs
      // for the second individual, etc).
      if (header[2] == '\x01') {
        // Get number of bytes.
        infile.seekg(0, infile.end);
        int num_bytes = infile.tellg();
        // Check if given dimensions match the file.
        if (n * p <= (num_bytes - 3) * 4) {
          // Iterate over rows indexes.
          for (int idx_i = 0; idx_i < size_i; idx_i++) {
            // Iterate over columns indexes.
            for (int idx_j = 0; idx_j < size_j; idx_j++) {
              // Reduce two-dimensional index to one-dimensional index with the mode.
              int which_pos = (j[idx_j] * n) + i[idx_i];
              // Every byte encodes 4 genotypes, find the one of interest.
              int which_byte = std::floor(which_pos / 4);
              // Find genotype in byte.
              int which_genotype = (which_pos % 4) * 2;
              // Read in the whole byte.
              infile.seekg(which_byte + 3);
              char *genotypes = new char[1];
              infile.read(genotypes, 1);
              // Remove the other genotypes by shifting the genotype of interest
              // to the end of the byte and masking with 00000011.
              int genotype = genotypes[0] >> which_genotype & 3;
              // Remap genotype value.
              int mapping;
              if (genotype == 0) {
                mapping = 0; // homozygous AA
              } else if (genotype == 3) {
                mapping = 2; // homozygous BB
              } else if (genotype == 1) {
                mapping = 1; // heterozygous AB
              } else if (genotype == 2) {
                mapping = NA_INTEGER; // missing
              } else {
                Rcpp::stop("Invalid genotype.");
              }
              out(idx_i, idx_j) = mapping;
            }
          }
        } else {
          Rcpp::stop("n or p does not match the dimensions of the file.");
        }
      } else {
        Rcpp::stop("Individual-major mode is not supported.");
      }
    } else {
      Rcpp::stop("File is not a binary PED file.");
    }
  } else {
    Rcpp::stop("File not found.");
  }
  return out;
}
