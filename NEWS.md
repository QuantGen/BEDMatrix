# BEDMatrix 2.0.2

- Fix that BED files larger than 2GB cannot be loaded in Windows.
- Do not produce `-Wunused-function` warnings in exported `BEDMatrix.h` header.


# BEDMatrix 2.0.1

- Fix extra `;` outside of `compute_num_bytes_per_variant` function in exported
  `BEDMatrix.h` header (found by CRAN incoming checks).


# BEDMatrix 2.0.0

- Extractions can be canceled with `Control-c`.
- Reimplement package in C, drop Rcpp.
- Genotypes can be extracted on the C level in other packages by adding
  `BEDMatrix` to the `LinkingTo` field of the `DESCRIPTION` file and including
  the `BEDMatrix.h` in your C code. The header file contains the `BEDMatrix`
  struct and the following functions: `compute_num_bytes_per_variant()`,
  `extract_genotype_linear()`, `extract_genotype_cartesian()`, and
  `recode_genotype()`.


# BEDMatrix 1.6.1

- Re-add C++11 requirement to avoid warnings from Boost headers on Windows.


# BEDMatrix 1.6.0

- Follow [Bioconductor S4 practices][2]. If you have used `new()` to create
  BEDMatrix instances, please use the constructor `BEDMatrix()` instead.
- Fix non-character rownames or colnames when data.table package is installed.
- Drop C++11 requirement.
- Update citation instructions.


# BEDMatrix 1.5.0

- Add `simple_names` option to BEDMatrix initialization function (contributed
  by Alex Ochoa).
- Fix single indexing for index values greater than `.Machine$integer.max`
  (previously generated `NA`s).
- Minor performance improvements.


# BEDMatrix 1.4.1

- Specify fileset name in messages during instantiation.


# BEDMatrix 1.4.0

- Support subsetting by `NA`.
- Use [crochet package][1] for subsetting
- Fix minor subsetting bugs
- Update example
- Remove intermediate S3 `BEDMatrix` type


# BEDMatrix 1.3.0

- Speed up character subsetting.
- Speed up detection of `n`, `p`, `rownames`, and `colnames` during
  initialization if `data.table` package is installed.
- Support `str` function.


# BEDMatrix 1.2.0

- Support `path` without extension (like PLINK).
- Add `path` attribute (to reattach instance after saving it to RData).
- Add `length` method.
- Add `as.matrix` method.
- Add `is.matrix` method.
- Store dimensions in S3 wrapper as `dims` attribute to allow for faster
  recreation when saved.
- Fix bug that modified `i` and `j` when subsetting.


# BEDMatrix 1.1.0

- Restore cross-platform compatibility by dropping the system dependency on
  `Boost.IOStreams` in favor of `Boost.Interprocess` which provides header-only
  memory-mapping and is therefore supported by the `BH` package.


# BEDMatrix 1.0.1

- Ensure that the same C compiler and compiler flags are used in the configure
  tests as when compiling R.
- Improve messages in configure script to distinguish between Boost headers and
  Boost libraries.


# BEDMatrix 1.0.0

Initial release.

[1]: https://CRAN.R-project.org/package=crochet
[2]: https://bioconductor.org/help/course-materials/2017/Zurich/S4-classes-and-methods.html
