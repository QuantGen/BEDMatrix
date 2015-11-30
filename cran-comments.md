## Release summary

Sorry for the frequent updates, but this will be the last one for a while. This
release drops the system dependency on `Boost.IOStreams` in favor of
`Boost.Interprocess` which provides header-only memory-mapping and is therefore
supported by the `BH` package. This should solve all the problems that were
encountered previously when installing the package.

---

## Test environments

* Local Arch Linux install, R 3.2.2
* Ubuntu 12.04 on Travis CI, R. 3.2.2
* win-builder, R-release and R-devel

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies

There are currently no downstream dependencies for this package.
