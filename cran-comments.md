## Release summary

As requested by Prof Ripley, I
* ensured that the same C compiler and compiler flags are used in the configure
  tests as when compiling R.
* improved the messages in the configure script to distinguish between Boost
  headers and Boost libraries.

Thanks for the feedback.

There are still two ERRORs on the package check page: the one on the
`r-release-osx-x86_64-mavericks` flavor is related to [Boost not being
installed](https://www.r-project.org/nosvn/R.check/r-release-osx-x86_64-mavericks/BEDMatrix-00install.html),
and the other one on the `r-devel-linux-x86_64-fedora-clang` flavor is maybe
related to using the wrong compiler, an issue that has been addressed in this
release.

---

This package is currently only supported on UNIX, as indicated by `OS_type:
unix` in the DESCRIPTION file.  The reason for this is that the package depends
on Boost.IOStreams, which is not header-only and therefore not available as
part of the BH package. I have not succeeded in setting up Boost on a Windows
machine yet. I have noted the external dependency in the `SystemRequirements`
field of the DESCRIPTION file. I also provided a configure script that tests
for the presence of Boost and Boost.IOStreams.

## Test environments

* Local Arch Linux install, R 3.2.2
* Ubuntu 12.04 on Travis CI, R. 3.2.2
* win-builder, R-release and R-devel

win-builder was used as a test environment, but both R-release and R-devel did
not "attempt to install this package on Windows" without any further ERRORS,
WARNINGs, or NOTEs.

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
