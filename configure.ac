AC_PREREQ([2.69])

AC_INIT(BEDMatrix, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))

AC_CONFIG_AUX_DIR([tools])

# Ensure same compiler and compiler flags as R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    echo "could not determine R_HOME"
    exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CXX=`"${R_HOME}/bin/R" CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)

m4_include([tools/m4/ax_boost_base.m4])
AX_BOOST_BASE([1.46.0],,AC_MSG_ERROR(Please install Boost.))

m4_include([tools/m4/ax_boost_iostreams.m4])
AX_BOOST_IOSTREAMS

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
