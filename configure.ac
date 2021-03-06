# Process this file with autoconf to produce a configure script.

# Initialization
AC_PREREQ([2.63])
AC_INIT([LatAnalyze],[3.2.1],[antonin.portelli@me.com],[LatAnalyze])
AC_CONFIG_AUX_DIR([.buildutils])
AC_CONFIG_SRCDIR([lib/Global.cpp])
AC_CONFIG_SRCDIR([utils/sample_read.cpp])
AC_CONFIG_SRCDIR([examples/exMathInterpreter.cpp])
AC_CONFIG_MACRO_DIR([.buildutils/m4])
AM_INIT_AUTOMAKE([1.11 -Wall -Werror foreign])
AM_SILENT_RULES([yes])
AC_CONFIG_HEADERS([config.h])
AM_CONDITIONAL([HAVE_AM_MINOR_LE_11],
               [test `automake --version | grep automake | awk -F '.' '{print $2}'` -le 11])
# Checks for programs
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LEX
AC_PROG_YACC
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])

LT_INIT

# Options
AC_ARG_WITH([gsl],
    [AS_HELP_STRING([--with-gsl=prefix],
		[try this for a non-standard install prefix of the GSL library])],
	[AM_CXXFLAGS="$AM_CXXFLAGS -I$with_gsl/include"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_gsl/lib"])
AC_ARG_WITH([minuit],
    [AS_HELP_STRING([--with-minuit=prefix],
		[try this for a non-standard install prefix of the Minuit2 library])],
    [AM_CXXFLAGS="$AM_CXXFLAGS -I$with_minuit/include"]
	[AM_LDFLAGS="$AM_LDFLAGS -L$with_minuit/lib"])
AC_ARG_WITH([nlopt],
    [AS_HELP_STRING([--with-nlopt=prefix],
        [try this for a non-standard install prefix of the NLopt library])],
    [AM_CXXFLAGS="$AM_CXXFLAGS -I$with_nlopt/include"]
    [AM_LDFLAGS="$AM_LDFLAGS -L$with_nlopt/lib"])
AC_ARG_WITH([hdf5],
    [AS_HELP_STRING([--with-hdf5=prefix],
    [try this for a non-standard install prefix of the HDF5 library])],
    [AM_CXXFLAGS="$AM_CXXFLAGS -I$with_hdf5/include"]
    [AM_LDFLAGS="$AM_LDFLAGS -L$with_hdf5/lib"])
AC_ARG_WITH([latcore],
    [AS_HELP_STRING([--with-latcore=prefix],
    [use this option for a non-standard install prefix of the LatCore library])],
    [AM_CXXFLAGS="$AM_CXXFLAGS -I$with_latcore/include"]
    [AM_LDFLAGS="$AM_LDFLAGS -L$with_latcore/lib"])
CFLAGS="$AM_CFLAGS $CFLAGS"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"

# Get compilers informations
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([C_COMP_VENDOR],["$ax_cv_c_compiler_vendor"],
			[vendor of C compiler that will compile the code])
AM_CONDITIONAL([CC_GNU],[test $ax_cv_c_compiler_vendor = "gnu"])
AM_CONDITIONAL([CC_INTEL],[test $ax_cv_c_compiler_vendor = "intel"])
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])
AX_COMPILER_VENDOR
AC_DEFINE_UNQUOTED([CXX_COMP_VENDOR],["$ax_cv_cxx_compiler_vendor"],
			[vendor of C++ compiler that will compile the code])
AM_CONDITIONAL([CXX_GNU],[test $ax_cv_cxx_compiler_vendor = "gnu"])
AM_CONDITIONAL([CXX_INTEL],[test $ax_cv_cxx_compiler_vendor = "intel"])
AX_GCC_VERSION
AC_DEFINE_UNQUOTED([GCC_VERSION],["$GCC_VERSION"],
			[version of gcc that will compile the code])
AX_GXX_VERSION
AC_DEFINE_UNQUOTED([GXX_VERSION],["$GXX_VERSION"],
			[version of g++ that will compile the code])

# Checks for libraries.
AC_CHECK_LIB([m],[cos],[],[AC_MSG_ERROR([libm library not found])])
AC_CHECK_LIB([gslcblas],[cblas_dgemm],[],
             [AC_MSG_ERROR([GSL CBLAS library not found])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],[AC_MSG_ERROR([GSL library not found])])
AC_CHECK_LIB([nlopt_cxx],[nlopt_create],
    [AC_DEFINE([HAVE_NLOPT],
    [1],
    [Define to 1 if you have the `NLopt' library (-lnlopt_cxx).])]
    [have_nlopt=true]
    [LIBS="$LIBS -lnlopt_cxx"],[])
AM_CONDITIONAL([HAVE_NLOPT], [test x$have_nlopt = xtrue])
AC_CHECK_LIB([hdf5_cpp],[H5Fopen],
    [AC_DEFINE([HAVE_HDF5],
    [1],
    [Define to 1 if you have the `HDF5' library (-lhdf5_cpp).])]
    [have_hdf5=true]
    [LIBS="$LIBS -lhdf5_cpp -lhdf5"],[],[-lhdf5])
AM_CONDITIONAL([HAVE_HDF5], [test x$have_hdf5 = xtrue])
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lMinuit2"
AC_MSG_CHECKING([for ROOT::Minuit2::BasicMinimumError in -lMinuit2]);
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM([#include <Minuit2/BasicMinimumError.h>],
	[ROOT::Minuit2::BasicMinimumError dummy(0)])],
	[LIBS="$LIBS -lMinuit2"]
	[AC_DEFINE([HAVE_MINUIT2],
				[1],
				[Define to 1 if you have the `Minuit2' library (-lMinuit2).])]
    [have_minuit=true]
	[AC_MSG_RESULT([yes])],
	[have_minuit=false]
    [AC_MSG_RESULT([no])])
AM_CONDITIONAL([HAVE_MINUIT], [test x$have_minuit = xtrue])
LDFLAGS=$SAVED_LDFLAGS
SAVED_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS -lLatCore"
AC_MSG_CHECKING([for LatCore::XmlReader in -lLatCore]);
AC_LINK_IFELSE(
	[AC_LANG_PROGRAM([#include <LatCore/XmlReader.hpp>],
    [LatCore::XmlReader dummy()])],
	[LIBS="$LIBS -lLatCore"]
	[AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([LatCore library not found])])
LDFLAGS=$SAVED_LDFLAGS

# Checks for header files.
AC_HEADER_STDC

AC_SUBST([LIBS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile lib/Makefile utils/Makefile physics/Makefile 
                 examples/Makefile])
AC_OUTPUT

echo "*********************************************"
echo "* ${PACKAGE_NAME} v${VERSION}" build options
echo "*********************************************"
echo "* C++ compiler     : $CXX"
echo "* HDF5 support     : `test x$HAVE_HDF5_TRUE = 'x'   && echo yes \
                                                          || echo no`"
echo "* Minuit minimizers: `test x$HAVE_MINUIT_TRUE = 'x' && echo yes \
                                                          || echo no`"
echo "* NLopt minimizers : `test x$HAVE_NLOPT_TRUE = 'x'  && echo yes \
                                                          || echo no`"
echo "*********************************************"



