# BEDMatrix 1.4.1

* Specify fileset name in messages during instantiation.


# BEDMatrix 1.4.0

* Support subsetting by `NA`.
* Use [crochet package](https://CRAN.R-project.org/package=crochet) for subsetting
* Fix minor subsetting bugs
* Update example
* Remove intermediate S3 `BEDMatrix` type


# BEDMatrix 1.3.0

* Speed up character subsetting.
* Speed up detection of `n`, `p`, `rownames`, and `colnames` during
  initialization if `data.table` package is installed.
* Support `str` function.


# BEDMatrix 1.2.0

* Support `path` without extension (like PLINK).
* Add `path` attribute (to reattach instance after saving it to RData).
* Add `length` method.
* Add `as.matrix` method.
* Add `is.matrix` method.
* Store dimensions in S3 wrapper as `dims` attribute to allow for faster
  recreation when saved.
* Fix bug that modified `i` and `j` when subsetting.


# BEDMatrix 1.1.0

* Restore cross-platform compatibility by dropping the system dependency on
  `Boost.IOStreams` in favor of `Boost.Interprocess` which provides header-only
  memory-mapping and is therefore supported by the `BH` package.


# BEDMatrix 1.0.1

* Ensure that the same C compiler and compiler flags are used in the configure
  tests as when compiling R.
* Improve messages in configure script to distinguish between Boost headers and
  Boost libraries.


# BEDMatrix 1.0.0

Initial release.
