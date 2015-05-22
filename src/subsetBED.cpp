#include <Rcpp.h>

#include <string>
#include <iostream>
#include <fstream>

//' Subsets a BED file.
//'
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix subsetBED(Rcpp::String path, int n, int p, Rcpp::IntegerVector i, Rcpp::IntegerVector j) {
  // Check if indexes are out of bounds
  if (Rcpp::is_true(Rcpp::any(i > n)) || Rcpp::is_true(Rcpp::any(j > p))) {
    Rcpp::stop("Invalid dimensions.");
  }
  // Keep sizes of i and j
  int size_i = i.size();
  int size_j = j.size();
  // Reserve output matrix
  Rcpp::IntegerMatrix out (size_i, size_j);
  // Open BED file
  std::ifstream infile (path.get_cstring(), std::ios::binary);
  if (infile) {
    char *header = new char[3];
    infile.read(header, 3);
    // Check magic number
    if (header[0] == '\x6C' && header[1] == '\x1B') {
      // Check individual-major mode
      if (header[2] == '\x01') {
        // Iterate over rows
        for (int idx_i = 0; idx_i < size_i; idx_i++) {
          int cur_i = i[idx_i];
          // Iterate over columns
          for (int idx_j = 0; idx_j < size_j; idx_j++) {
            int cur_j = j[idx_j];
            // Reduce two-dimensional index to one-dimensional index
            int whichPos = ((cur_i - 1) * p) + (cur_j - 1);
            // Every byte encodes 4 genotypes, find the one of interest
            int whichByte = std::floor(whichPos / 4);
            // Find genotype in byte
            int whichGenotype = (whichPos % 4) * 2;
            // Read in the whole byte
            infile.seekg(whichByte + 3);
            char *genotypes = new char[1];
            infile.read(genotypes, 1);
            // Remove the other genotypes
            int genotype = genotypes[0] >> whichGenotype & 3;
            // Remap genotype value
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
        Rcpp::stop("Individual major mode is not supported.");
      }
    } else {
      Rcpp::stop("File is not a binary PED file.");
    }
    infile.close();
  } else {
    Rcpp::stop("File not found.");
  }
  return out;
}
