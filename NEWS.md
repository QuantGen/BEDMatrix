# BEDMatrix 1.1.0.9000

* Add `length` method.
* Add `as.matrix` method.

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
