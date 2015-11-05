## Resubmission

This is the second resubmission.

I received feedback asking whether PED is an acronym. It is not, but a well
established file format in the genetics community, especially in human
genetics.

In this version I have:
- Added "(PLINK)" in the title to give better context.
- Added a link to the PLINK website in the Description.

In the previous resubmission I have:
- Fixed the mis-use of the LICENSE file by removing the MIT license text and
  replacing it with `YEAR` and `COPYRIGHT` fields.
- Changed the Description field to not start with "This package".
- Changed the Description field to explain that a PED file is one of the file
  formats of PLINK, a popular toolset for whole genome association analysis.

Thanks for the feedback.

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
