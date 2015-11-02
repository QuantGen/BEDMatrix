AC_PREREQ([2.69])

AC_INIT(BEDMatrix, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))

AC_CONFIG_AUX_DIR([tools])

AC_REQUIRE_CPP

m4_include([tools/m4/ax_boost_base.m4])
AX_BOOST_BASE([1.46.0],,AC_MSG_ERROR(Please install Boost.))

m4_include([tools/m4/ax_boost_iostreams.m4])
AX_BOOST_IOSTREAMS
if test "x$BOOST_IOSTREAMS_LIB" = "x"; then
    AC_MSG_ERROR(Please install the Boost.IOStreams library.)
fi

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT
