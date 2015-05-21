#include <Rcpp.h>

#include <string>
#include <iostream>
#include <fstream>

//' Subsets a BED file.
//'
//' @export
// [[Rcpp::export]]
int subsetBED(Rcpp::String path, int n, int p, int i, int j) {
  int out;
  if (i > n || j > p) {
    Rcpp::stop("Invalid dimensions.");
  }
  std::ifstream infile (path.get_cstring(), std::ios::binary);
  if (infile) {
    char *header = new char[3];
    infile.read(header, 3);
    // Check magic number
    if (header[0] == '\x6C' && header[1] == '\x1B') {
      // Check individual-major mode
      if (header[2] == '\x01') {
        // Reduce two-dimensional index to one-dimensional index
        int whichPos = ((i - 1) * p) + (j - 1);
        // Every byte encodes 4 genotypes, find the one of interest
        int whichByte = std::floor(whichPos / 4);
        // Find genotype in byte
        int whichGenotype = (whichPos % 4) * 2;
        // Read in the whole byte
        infile.seekg(whichByte, std::ios::cur);
        char *genotypes = new char[1];
        infile.read(genotypes, 1);
        // Remove the other genotypes
        int genotype = genotypes[0] >> whichGenotype & 3;
        // Remap genotype value
        if (genotype == 0) {
          out = 0; // homozygous AA
        } else if (genotype == 3) {
          out = 2; // homozygous BB
        } else if (genotype == 1) {
          out = 1; // heterozygous AB
        } else if (genotype == 2) {
          out = NA_INTEGER; // missing
        } else {
          Rcpp::stop("Invalid genotype.");
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
