#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Few bits where taken from pyconf, bearing this notice:
dnl Copyright 2012, 2013 Brandon Invergo <brandon@invergo.net>
dnl
dnl Copying and distribution of this file, with or without modification,
dnl are permitted in any medium without royalty provided the copyright
dnl notice and this notice are preserved.  This file is offered as-is,
dnl without any warranty.


AC_PREREQ([2.71])
AC_INIT([Gyoto],[2.0.0],[gyoto@sympa.obspm.fr])
gyoto_test_CXXFLAGS=${CXXFLAGS+set}
gyoto_test_CPPFLAGS=${CPPFLAGS+set}
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([include/GyotoKerrKS.h])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
LT_INIT([dlopen disable-static])
AC_CONFIG_HEADERS([config.h include/GyotoConfig.h])
AC_SUBST([VERBATIM], [])


# http://www.gnu.org/software/libtool/manual/html_node/Updating-version-info.html
# C: current; R: revision; A: age
# Linux library name is: libgyoto.so.(C-A).A.R
VERSINFO_C=9
VERSINFO_R=0
VERSINFO_A=0
VERSINFO=${VERSINFO_C}:${VERSINFO_R}:${VERSINFO_A}
AC_SUBST([VERSINFO])

## Needed for gcc-5, should be harmless for other compilers:
#
AX_APPEND_COMPILE_FLAGS([-P], [CPPFLAGS])

## Ensure local build does not have a release soname by default
#
FEATURES=""
AC_ARG_ENABLE([release],
  [AS_HELP_STRING([--enable-release],
     [enable a release build (remove -unreleased from SONAME)])])
AS_IF([test "x$enable_release" != "xyes"],
  [FEATURES+="-unreleased"])

## An additional directory to look for plug-ins
#  Default /usr/local/lib/gyoto, only if ${prefix} is not /usr/local
AC_ARG_ENABLE([localpkglibdir],
  [AS_HELP_STRING([--enable-localpkglibdir[[=PKGLIBDIR]]],
     [additional directory to look for plug-ins [/usr/local/lib/gyoto, only if PREFIX is not /usr/local]. If specified, Gyoto will look for additional plug-ins in LOCALPKGLIBDIR/GYOTO_SOVERS and  LOCALPKGLIBDIR before looking in PKGLIBDIR/GYOTO_SOVERS and PKGLIBDIR])],
  [],
  [AS_IF([test x${prefix} == xNONE || test "x${prefix}" == x/usr/local],
     [enable_localpkglibdir=no],
     [enable_localpkglibdir=yes])])
AS_IF([test x$enable_localpkglibdir != "xno"],
  [AS_IF([test x$enable_localpkglibdir == "xyes"],
     [localpkglibdir="/usr/local/lib/"${PACKAGE}],
     [localpkglibdir=$enable_localpkglibdir]
  )],
  [localpkglibdir=""])
AM_CONDITIONAL(HAVE_LOCALPKGLIBDIR, [test "x$localpkglibdir" != "x"])
AC_SUBST(localpkglibdir)

## SIGFPE delivery
#
AC_ARG_WITH([fenv],
  [AS_HELP_STRING([--without-fenv],
     [do not use fenv.h (see BUGS.md)],)],
   [],
   [with_fenv=auto])
AS_IF([test "x$with_fenv" != xno],
  [AC_MSG_CHECKING([for fenv.h])
   AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <fenv.h>]], [[feenableexcept(FE_DIVBYZERO);]])],[AC_DEFINE([HAVE_FENV_H], [1], [Define to 1 to use fenv.h.])
      HAVE_FENV_H=1
      AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
  ])
AS_IF([test "x$with_fenv" = xyes && test "x$HAVE_FENV_H" != x1],
  [AC_MSG_ERROR(fenv.h required but not usable")])

## IS SIZE_T EQUIVALENT TO UNSIGNED LONG?
#
AC_MSG_CHECKING(whether size_t is the same as unsigned long)

AC_ARG_WITH([size_t],
  [AS_HELP_STRING([--with-size_t=same|distinct|auto],
     [consider size_t as the same as or distinct from unsigned long (auto). Should normally be specified only if cross-compiling.])],
  [AS_IF([test "x$with_size__t" != xsame && test "x$with_size__t" != xdistinct && test "x$with_size__t" != xauto],
     [AC_MSG_ERROR([with-size_t should be one of same, distinct or auto])])
  ],
  [with_size__t=auto]
)

AS_IF([test "x$with_size__t" = xauto],
  [AC_RUN_IFELSE([AC_LANG_SOURCE([[
      #include <sys/types.h>
      typedef void (* set_size_t_t) (size_t val);
      typedef void (* set_unsigned_long_t) (unsigned long);
      size_t size_t_value;
      void set_size_t (size_t v) {size_t_value = v;}
        int main(int, char**) {
        set_unsigned_long_t ulsetter = &set_size_t;
        (*ulsetter)(25);
        if (size_t_value != 25) return 1;
        return 0;
      }
     ]])],[with_size__t="same"],[with_size__t="distinct"],[])
  ]
)

AS_IF([test "x$with_size__t" = xsame],
  [AC_DEFINE([GYOTO_SIZE__T_IS_UNSIGNED_LONG], [1], [Define to 1 if size_t and unsigned long are the same])]
)

AC_MSG_RESULT($with_size__t)

# PYTHON STUFF
AC_ARG_WITH([python],
	[AS_HELP_STRING([--with-python@<:@=interpreter@:>@],[force (or disable) python support, possibly specifying interpreter])],
	[
	 case $with_python in
	   yes | no) ;;
	   *)
	       PYTHON=`which $with_python` || AC_MSG_ERROR([$with_python: command not found])
	       with_python=yes
	 esac
	],
	[])

PC_PROG_PYTHON([python3], [3.0], [4.0])
PC_PYTHON_PROG_PYTHON_CONFIG
PC_PYTHON_CHECK_HEADERS(
 [
  PC_PYTHON_CHECK_INCLUDES
  PC_PYTHON_CHECK_CFLAGS
  PC_PYTHON_CHECK_LDFLAGS
  PC_PYTHON_CHECK_LIBS
  PC_PYTHON_CHECK_SITE_DIR
  PC_PYTHON_SITE_PACKAGE_DIR
  PC_PYTHON_CHECK_EXEC_DIR
  PC_PYTHON_CHECK_EXTENSION_SUFFIX
  PC_PYTHON_CHECK_EXEC_DIR
],[PYTHON=""])
AS_IF([test "x$with_python" = xno], [PYTHON=""])
AS_IF([test "x$with_python" = xyes && test "x$PYTHON" = x],
  [AC_MSG_ERROR([Python requested but no supported version found (>=3.0, <4.0)])])

AC_ARG_WITH([virtualenv],
	[AS_HELP_STRING([--with-virtualenv], [install to a Python virtualenv])],
	[],
	[with_virtualenv=no])
AS_IF([test "x$with_virtualenv" != xno],
	    [AC_CHECK_PROGS([VIRTUALENV], [virtualenv virtualenv3 virtualenv2], [no])]
	    [AS_IF([test "x$VIRTUALENV" = xno],
	    [AC_MSG_FAILURE(
		[--with-virtualenv given but virtualenv could not be found])])],
	[AC_SUBST([VIRTUALENV], [no])])
AC_ARG_VAR([VIRTUALENV_FLAGS], [flags to pass to the virtualenv command])

AX_PKG_SWIG([2.0])

# DONE WITH PYTHON STUFF

# Checks for programs.

############### CHECK FOR OPTIONAL FEATURES #########################
## ADVANCED INTEGRATORS
# Check for C++11
AC_ARG_ENABLE([c++11],
  [AS_HELP_STRING([--disable-c++11],
     [disable C++11 features (including advanced integrators)])])
AS_IF([test "x$enable_c__11" != "xno"],
  [AX_CXX_COMPILE_STDCXX_11(, [optional])]
)
AS_IF([test "x$HAVE_CXX11" != "x1" && test "x$enable_c__11" = "xyes"],
  [AC_MSG_ERROR([C++11 requested but not found])])

AC_ARG_WITH([mpi], [AS_HELP_STRING([--with-mpi],
    [compile with MPI (parallelization) support. If none is found,
    MPI is not used. Default: auto])
],,[with_mpi=auto])

# Check for boost

have_boost=no
AS_IF([test "x$with_boost" != "xno"],
  [have_boost=yes
   BOOST_REQUIRE([1.53.1], [have_boost=no])
   AS_IF([test "x$have_boost" == "xyes"],
     [BOOST_FIND_HEADER([boost/multiprecision/cpp_dec_float.hpp],
        [AC_MSG_WARN([Boost.multiprecision not found])
	 have_boost=partial],
        [])
      BOOST_FIND_HEADER([boost/array.hpp],
        [AC_MSG_WARN([boost/array.hpp not found])
	 have_boost=partial],
        [])
      BOOST_FIND_HEADER([boost/numeric/odeint/stepper/controlled_step_result.hpp],
        [AC_MSG_WARN([Boost.odeint not found])
	 have_boost=partial],
        [])
     ])
   ],
   [have_boost=no]
)

AS_IF([test "x$HAVE_CXX11" == "x1" &&  test "x$ac_cv_header_boost_numeric_odeint_stepper_controlled_step_result_hpp" == xyes],
  [AC_DEFINE([GYOTO_HAVE_BOOST_INTEGRATORS], [1], [Define to 1 if C++11 and Boost.odeint are both available]
  )]
)

AS_IF([test "x$have_boost" != "xyes"],
  [FEATURES+="-noboost"]
)

# If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.

have_mpi=no
AC_CHECK_PROGS([MPICXX], [mpic++ mpicxx])
AS_IF([test "x$MPICXX" != x && test "x$have_boost" == xyes && test "x$with_mpi" != xno],
  [
   ax_save_CXX=$CXX
   ax_save_CPP=$CXXCPP
   CXX=$MPICXX
   CXXCPP="$MPICXX -E"
   BOOST_FIND_HEADER([boost/mpi/environment.hpp],
     [have_mpi=no
      CXX=$ax_save_CXX
      AS_IF([test "x$with_mpi" == xyes],
        [ AC_MSG_ERROR([MPI features require boost.mpi])
      ])],
     [
      BOOST_FIND_LIBS([mpi], [mpi], [], [boost/mpi/environment.hpp],
        [boost::mpi::environment()], [])
      BOOST_FIND_LIBS([serialization], [serialization], [],
        [boost/serialization/base_object.hpp],
        [], [])
      have_mpi=yes
    ])
  ]
)

AS_IF([test "x$have_mpi" == xyes],
[ AC_SUBST(MPICXX)
  AC_DEFINE([HAVE_MPI], [1], [Define to 1 if you have MPI])],
[ AS_IF([test "x$with_mpi" == xyes],
  [ AC_MSG_ERROR(MPI requested but not found)
  ],
  [FEATURES+="-nompi"
   MPICXX=""])
])

AM_CONDITIONAL(HAVE_MPI, [test "x$have_mpi" == xyes])

AC_PROG_MKDIR_P
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
PKG_PROG_PKG_CONFIG([])
AC_CHECK_PROGS([DOXYGEN], [doxygen])
AC_CHECK_PROGS([PDFLATEX], [pdflatex])
AC_CHECK_PROGS([BIBTEX], [bibtex])
AC_CHECK_PROGS([KPSEWHICH], [kpsewhich])

# Provide additional defaults for user variables

AX_APPEND_COMPILE_FLAGS([-rdynamic], [LDFLAGS])

AC_ARG_ENABLE([native],
  [AS_HELP_STRING([--enable-native],[append "-march=native" to CXXFLAGS if supported])],
  [AS_IF([test "x$enable_native" = "xyes"],
     [AX_APPEND_COMPILE_FLAGS([[-march=native]])])])

AC_ARG_ENABLE([hardening],
  [AS_HELP_STRING([--enable-hardening],
     [add warnings and hardening flags if supported])],
  [AS_IF([test "x$enable_hardening" = "xyes"],
     [AX_APPEND_COMPILE_FLAGS([-pedantic -Wall -W -Wundef -Wshadow -Wcast-qual \
        -Wcast-align -Wconversion -Winline -Wabi -Woverloaded-virtual \
        -Wno-long-long -Wformat -Wformat-security -Werror=format-security \
        -fstack-protector "--param ssp-buffer-size=4"])
      AX_APPEND_COMPILE_FLAGS([-D_FORTIFY_SOURCE=2], [CPPFLAGS])
      AX_APPEND_COMPILE_FLAGS([-fPIE], [CXXFLAGS_EXEC])
      AX_APPEND_LINK_FLAGS([-pie], [LDFLAGS_EXEC], [-fPIE])
     ])
  ]
)

AC_SUBST([CXXFLAGS_EXEC])
AC_SUBST([LDFLAGS_EXEC])

AC_ARG_ENABLE([deprecated],
  [AS_HELP_STRING([--enable-deprecated],
     [enable deprecated method names (needing this is a bug)])],
  [AS_IF([test "x$enable_deprecated" = "xno"],
     [AX_APPEND_COMPILE_FLAGS([-DGYOTO_NO_DEPRECATED], [CPPFLAGS])])],
  [AX_APPEND_COMPILE_FLAGS([-DGYOTO_NO_DEPRECATED], [CPPFLAGS])]
)


# Guess *LD_LIBRARY_PATH variable for make check
AC_SUBST([target_os])
SYS=other
DYLIB_VAR=LD_LIBRARY_PATH
case "x$target_os" in
     xdarwin*|xmacos)
	DYLIB_VAR=DYLD_LIBRARY_PATH
	SYS=Darwin
	;;
     *) ;;
esac
plugin_sfx=so
AC_DEFINE_UNQUOTED([GYOTO_PLUGIN_SFX],["${plugin_sfx}"],
  [Suffix for Gyoto plug-ins])
AC_SUBST([GYOTO_PLUGIN_SFX],[${plugin_sfx}])
AC_SUBST([DYLIB_VAR])
# Used in yorick/stdplug/Makefile
AC_SUBST([SYS])
AC_SUBST([plugin_sfx])

# Checks for libraries.
AC_CHECK_LIB([dl], [dlopen])

# Checks for header files.
# We are not using those tests to far
# AC_CHECK_HEADERS([float.h stddef.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
# other functions we may need to check for according to autoscan:
# floor getcwd pow sqrt strdup
AC_CHECK_FUNCS([sincos])

## YORICK
# Check for Yorick
yorick_exe=yorick
YORICK=
AC_MSG_CHECKING([whether to use yorick])
AC_ARG_WITH([yorick],
  [AS_HELP_STRING([--with-yorick[[=yorick-executable]]],[force using yorick and optionnally specify which one [default=yorick,
      only if found]])],
  [AS_IF([test "x$with_yorick" == "xno"],
     [AC_MSG_RESULT([no])],
     [AS_IF([test "x$with_yorick" != "xyes"],
        [YORICK=$with_yorick],
        [YORICK=$yorick_exe])
      AS_IF([$YORICK -batch],
        [AC_MSG_RESULT([$YORICK])],
        [AC_MSG_ERROR([$YORICK does not work])])
     ])
   AC_SUBST([YORICK])
  ],
  [AC_CHECK_PROGS([YORICK],$yorick_exe)]
)

AM_CONDITIONAL(HAVE_YORICK, [test "x$YORICK" != "x"])

# Only if YORICK available, check for install directories
AS_IF([test "x$YORICK" != "x"],
  [AC_MSG_CHECKING([install path for Yorick platform-dependent files])
   Y_EXE_HOME=`echo write,format=\"%s\",Y_HOME\; quit | ${YORICK}| tail -1`
   Y_EXE_SITE=`echo write,format=\"%s\",Y_SITE\; quit | ${YORICK}| tail -1`
   Y_INST_HOME=
   Y_INST_SITE=
   AC_ARG_WITH([yorick-home],
     [AS_HELP_STRING([--with-yorick-home=DIR],[architecture-dependent Yorick files [default: the
         yorick directory if it is under ${exec_prefix}, else
         ${exec_prefix}/lib/yorick]])],
     [Y_INST_HOME=${with_yorick_home}
      Y_INST_SITE=${with_yorick_home}
     ],
     [AS_IF([echo ${Y_EXE_HOME} | grep ${exec_prefix}],
        [Y_INST_HOME=Y_EXE_HOME
	 Y_INST_SITE=Y_EXE_SITE],
	[Y_INST_HOME='${exec_prefix}/lib/yorick'
	 Y_INST_SITE='${prefix}/lib/yorick']
      )])
   AC_MSG_RESULT([$Y_INST_HOME])
   AC_SUBST([Y_INST_HOME])

   AC_MSG_CHECKING([install path for Yorick platform-independent files])
   AC_ARG_WITH([yorick-site],
     [AS_HELP_STRING([--with-yorick-site=DIR],[architecture-independent Yorick files
	 [default: the yorick directory if it is under ${prefix}, else
	  ${prefix}/lib/yorick]])
     ],
     [Y_INST_SITE=${with_yorick_site}])
   AC_MSG_RESULT([$Y_INST_SITE])
   AC_SUBST([Y_INST_SITE])

  ]
)

## DEBUGGING INFORMATION
# Enable debugging by default
GYOTO_DEBUG_ENABLED=yes
AC_ARG_ENABLE([debugging],
  [AS_HELP_STRING([--disable-debugging],
     [remove debugging code for hypothetical speed gain])],
  [AS_IF([test "x$enable_debugging" == "xno"],
     [GYOTO_DEBUG_ENABLED=0])
  ])
AS_IF([test "x$GYOTO_DEBUG_ENABLED" == xyes],
  [AC_DEFINE([GYOTO_DEBUG_ENABLED],
     [1],
     [Enable debugging information gathering code for putative speed gain])
  ]
)

## POSIX THREADS
# Check for pthreads
use_pthreads=yes
AC_ARG_ENABLE([pthreads],
  [AS_HELP_STRING([--disable-pthreads],
     [disable POSIX threads parallelization])],
  [AS_IF([test "x$enable_pthreads" = "xno"], [use_pthreads=no])])
AS_IF([test "x$use_pthreads" == "xyes"], [AX_PTHREAD],
      [FEATURES+="-nopthreads"])

## UDUNITS
# Check for UDUNITS
pkg_requires=
pkg_cflags=
pkg_libs=

AC_ARG_WITH([udunits-inc],
  [AS_HELP_STRING([(--with-udunits-inc=)],
     [deprecated, use --with-udunits-headers])],
  [with_udunits_headers=$withval
   AC_MSG_WARN([--with-udunits-inc option is deprecated, use --with-udunits-headers instead])
  ]
)

AC_ARG_WITH([udunits-lib],
  [AS_HELP_STRING([(--with-udunits-lib=)],
     [deprecated, use --with-udunits-libs])],
  [with_udunits_libs=${withval}
   AC_MSG_WARN([--with-udunits-lib option is deprecated, use --with-udunits-libs instead])
  ]
)


GYOTO_ARG_LIB(udunits, udunits, -ludunits2,
	      [#include <udunits2.h>],
	      [ut_set_error_message_handler(ut_ignore);])
AS_IF([test "x$have_udunits" == "xno"],
  [FEATURES+="-nounits"]
)

## XERCES
# Check for Xerces-c
GYOTO_ARG_LIB(xerces, xerces-c, -lxerces-c,
	      [#include <xercesc/util/PlatformUtils.hpp>],
	      [xercesc::XMLPlatformUtils::Initialize();])
AS_IF([test "x$have_xerces" == "xno"],
  [FEATURES+="-noxml"
   AS_IF([test "x$have_mpi" == "xyes"],
     [AC_MSG_ERROR([The MPI feature set requires Xerces. Please reconfigure with --without-mpi or with --with-xerces.])])
  ]
)

## EIGEN
# Check for EIGEN
GYOTO_ARG_LIB(eigen, eigen3, ,
              [#include <Eigen/Dense>],
              [Eigen::MatrixXd m(2,2);])
AS_IF([test "x$have_ceigen" == "xno"],
  [FEATURES+="-noeigen"]
)

## CFITSIO
# Check for CFITSIO
GYOTO_ARG_LIB(cfitsio, cfitsio, -lcfitsio,
	      [#include <fitsio.h>
	       #include <string>],
	      [fitsfile* fptr=NULL; std::string file=""; int status=0;
	       fits_open_file(&fptr, file.c_str(), 0, &status);])
AS_IF([test "x$have_cfitsio" == "xno"],
  [FEATURES+="-nofits"]
)

## ARBLIB
# Check for ARBLIB
have_2f1=no
# If AEAE if explicitely requested, skip ARBLIB checks.
case "x$with_aeae" in
     x)    ;;
     xno)  ;;
     x*)  with_arblib=no ;;
esac
AC_MSG_CHECKING(whether arblib library name is provided)
AS_IF([test "x$with_arblib_ldflags" == "x"],
  [arblib_ldflags=-lflint-arb],
  [with_arblib=yes
   arblib_ldflags=$with_arblib_ldflags]
)
AC_MSG_RESULT([$arblib_ldflags])
GYOTO_ARG_LIB(arblib, arblib, $arblib_ldflags,
	      [#include <acb_hypgeom.h>],
	      [double hypergeom = 0., rad=0., kappaIndex=5., thetae=1;
  	 acb_t FF, aa, bb, cc, zed;
	 acb_init(FF);
	 acb_init(aa);
	 acb_init(bb);
	 acb_init(cc);
	 acb_init(zed);
	 acb_set_d_d(aa,   kappaIndex-1./3.,  0.);
	 acb_set_d_d(bb,   kappaIndex+1.,     0.);
	 acb_set_d_d(cc,   kappaIndex+2./3.,  0.);
	 acb_set_d_d(zed, -kappaIndex*thetae, 0.);
	 slong prec=53;
	 acb_hypgeom_2f1(FF, aa, bb, cc, zed, ACB_HYPGEOM_2F1_AC, prec);
	 hypergeom = arf_get_d(&acb_realref(FF)->mid, ARF_RND_NEAR);
	 rad = mag_get_d(&acb_realref(FF)->rad);
	 acb_clear(FF);
	 acb_clear(aa);
	 acb_clear(bb);
	 acb_clear(cc);
	 acb_clear(zed);])
AC_ARG_WITH([arblib-ldflags],
  [AS_HELP_STRING([--with-arblib-ldflags=-l<libname>],
     [linker flag for arblib (-lflint-arb on Debian derivatives, -larb on Arch Linux...) possibly with dependencies (e.g. "-larb -lflint").])]
)
AS_IF([test "x$have_arblib" == "xyes"],
  [have_2f1=yes
   with_aeae=no
  ],
  [arblib_ldflags=
  ]
)

## AEAE
# Check for AEAE
AC_MSG_CHECKING(whether aeae should be used)
aeae_headers=
AC_ARG_WITH(
  [aeae],
  [AS_HELP_STRING( [--with-aeae[[=/path/to/AEAE]]],
     [Force using AEAE, optionally specifying full path to location of source code, or force disable it.])],
  [
   AS_IF([test "x$with_aeae" == "xno"],
     [AC_MSG_RESULT(no)],
     [test "x$with_aeae" != "xyes"],
     [aeae_headers=$with_aeae
      case "x$aeae_headers" in
	   x\-I*) ;;
	   x) ;;
	   x*) aeae_headers=-I$aeae_headers ;;
      esac
      with_aeae=yes
      AC_MSG_RESULT(yes)],
     [AC_MSG_RESULT(yes)])
  ],
  [with_aeae=check
   AC_MSG_RESULT(probe)
  ]
)
have_aeae=no
AS_IF([test "x$with_aeae" != "xno"],
  [
   AS_IF([test "x$have_aeae" == "xno"],
     [AC_MSG_CHECKING([for aeae])
      TMPCPPFLAGS=$CPPFLAGS
      CPPFLAGS="$TMPCPPFLAGS $aeae_headers"
      AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include <complex>
	 #include <iostream>
	 #define SIGN(a) (((a) < 0) ? (-1) : (1))
	 using namespace std;
	 #include "complex_functions.H"
	 #include "hyp_2F1.cpp"]], [[double hypergeom = 0., kappaIndex=5., thetae=1;
	 complex<double> aa=kappaIndex-1./3., bb=kappaIndex+1.,
	 cc=kappaIndex+2./3., zed=-kappaIndex*thetae;
	 const complex<double> FF = hyp_2F1(aa,bb,cc,zed);
	 hypergeom = FF.real();
	]])],[have_aeae=yes
         AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no)])
      CPPFLAGS=$TMPCPPFLAGS
      AEAE_CFLAGS="$aeae_headers"
      pkg_cflags="${pkg_cflags} ${AEAE_CFLAGS}"
     ]
   )
  ],
  []
)
AS_IF([test "x$have_aeae" == "xno"],
  [AS_IF([test "x$with_aeae" == "xyes"],
     [AC_MSG_ERROR([aeae requested but not found])
     ],
     []
   )
  ],
  [have_2f1=yes
   AC_DEFINE([GYOTO_USE_AEAE],
     [1],
     [Define to 1 if you have aeae])
  ]
)
AM_CONDITIONAL([HAVE_AEAE], [test "x$have_aeae" == "xyes"])
AM_CONDITIONAL([HAVE_2F1], [test "x$have_2f1" == "xyes"])
AC_SUBST([AEAE_CFLAGS])


## LORENE
# Check for LORENE
AC_MSG_CHECKING([for LORENE])
AC_ARG_VAR([MAKE_S],
  [make command which does not print directories [[make -w]]])
AS_IF([test "x$MAKE_S" = "x"],
  [MAKE_S="make -s"])
AC_ARG_VAR([HOME_LORENE],
  [path to Lorene build directory (see also --with-lorene)])
# first check whether LORENE is required and get HOME_LORENE
require_lorene=no
AC_ARG_WITH([lorene],
  [AS_HELP_STRING([--with-lorene[[=${HOME_LORENE}]]],
     [force using the optional LORENE library, possibly specifying HOME_LORENE
      [default=${HOME_LORENE}]])],
  [AS_IF([test "x$with_lorene" == "xno"],
     [HOME_LORENE=""
      AC_MSG_RESULT([no])
     ],
     [test "x$with_lorene" == "xyes"],
     [require_lorene=yes
      AS_IF([test "x$HOME_LORENE" != "x"],
      [AC_MSG_RESULT($HOME_LORENE)],
      [AC_MSG_ERROR([--with-lorene passed with no path: please specify LORENE
                    path either using the HOME_LORENE environment variable
                    or --with-lorene=<path>.])])
     ],
     [require_lorene=yes
      HOME_LORENE="$with_lorene"
      AC_MSG_RESULT($HOME_LORENE)
     ])
  ],
  [AS_IF([test "x$HOME_LORENE" != "x"],
     [AC_MSG_RESULT($HOME_LORENE)],
     [AC_MSG_RESULT([no])])
  ]
)
# then check whether LORENE is actually installed in HOME_LORENE
# finally output HOME_LORENE
AC_SUBST([HOME_LORENE])
AM_CONDITIONAL([HAVE_LORENE], [test "x$HOME_LORENE" != "x"])
AS_IF([test "x$HOME_LORENE" != "x"],
  [
   AS_IF([test ! -r ${HOME_LORENE}/local_settings],
     [AC_MSG_ERROR([${HOME_LORENE}/local_settings not found])])
cat >conftest.make <<_ACEOF
include \$(HOME_LORENE)/local_settings
conftest.make:
	touch conftest.make

%: conftest.make
	@echo \$(\$@)
_ACEOF
   LORENELDFLAGS="-L$HOME_LORENE/Lib -llorene"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_LAPACK | tr -d '\n\r'`"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_GSL | tr -d '\n\r'`"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_PGPLOT | tr -d '\n\r'`"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_LORENE | tr -d '\n\r'`"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_CXX | tr -d '\n\r'`"
   LORENECPPFLAGS=`$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make INC | tr -d '\n\r'`
   LORENECPPFLAGS="$LORENECPPFLAGS"
   rm conftest.make
   AC_SUBST([LORENECPPFLAGS])
   AC_SUBST([LORENELDFLAGS])
   pkg_cflags="${pkg_cflags} ${LORENECPPFLAGS}"
   pkg_libss="${pkg_libss} ${LORENELDFLAGS}"

  ]
)

## DOCUMENTATION
# disable doc
AC_ARG_ENABLE([doc],
  [AS_HELP_STRING([--disable-doc],
     [disable building documentation (automatically)])],
  [AS_IF([test "x$enable_doc" == "xyes"],
     [mkdoc=yes],
     [test "x$enable_doc" == "xno"],
     [mkdoc=no],
     [AC_MSG_ERROR([Please set --disable-doc or --enable-doc only to ``yes'' or ``no''])])
  ],
  [mkdoc=yes])

# check if doxygen is available
AM_CONDITIONAL([MKREF], [test "x$DOXYGEN" != "x"])

# check if all the necessary LaTeX tools are available
AM_CONDITIONAL([MKGUIDE], [test "x$KPSEWHICH" != "x" && ${KPSEWHICH} kluwer.bst && test "x$PDFLATEX" != "x" && test "x$BIBTEX" != "x"])

####################   FINALIZATION   ###################################
# Finalization: select output files, subdirs etc.
AC_SUBST([pkg_requires])
AC_SUBST([pkg_cflags])
AC_SUBST([pkg_libs])
mySUBDIRS="lib"
build_exec=no
AS_IF([test "x$have_mpi" == xyes], [build_exec=yes])
AS_IF([test "x$use_xerces" != "x" && test "x$use_cfitsio" != "x"], [build_exec=yes])
AS_IF([test "x$build_exec" == xyes], [mySUBDIRS="$mySUBDIRS bin"])

AM_CONDITIONAL([BUILD_GYOTO], [test "x$use_xerces" != "x" && test "x$use_cfitsio" != "x"])

sovers=$(( VERSINFO_C - VERSINFO_A ))
sovers+=$FEATURES
AC_SUBST([FEATURES])
AC_SUBST([sovers])
AC_DEFINE_UNQUOTED([GYOTO_SOVERS], ["${sovers}"], [Gyoto ABI version])

AS_IF([test "x$PYTHON" != "x"],
 [AC_CONFIG_SUBDIRS([plugins/python])
  mySUBDIRS="$mySUBDIRS python"
 ])

AC_CONFIG_FILES([Makefile
                 bin/Makefile
                 doc/Makefile
                 doc/doxyfile
                 doc/user_guide/Makefile
                 lib/Makefile
                 lib/gyoto.pc
                 lib/gyoto-uninstalled.pc
                 ])
AS_IF([test "x$YORICK" != "x"],
  [mySUBDIRS="$mySUBDIRS yorick"
   AC_CONFIG_FILES([yorick/Makefile],
 		   [oldpwd=`pwd`
		    cd yorick
		    $YORICK -batch make.i
		    cd $oldpwd],
		   YORICK=$YORICK)
   AC_CONFIG_FILES([yorick/stdplug/Makefile],
 		   [oldpwd=`pwd`
		    cd yorick/stdplug
		    $YORICK -batch make.i
		    cd $oldpwd],
		   YORICK=$YORICK)
   AC_CONFIG_FILES([yorick/gyoto.info
		    yorick/yorick1
		    yorick/setpaths.i
		   ],
		   [chmod a+x yorick/yorick1])
   AC_CONFIG_FILES([gyoto],
		   [chmod a+x gyoto])
  ]
)
AS_IF([test "x$PYTHON" != "x"],
  [
   AC_CONFIG_FILES([python/Makefile python/setup.py])
  ]
)
AS_IF([test "x$DOXYGEN" != "x" && test "x$mkdoc" == "xyes"],
  [mySUBDIRS="$mySUBDIRS doc"
  ]
)
AC_SUBST([mySUBDIRS])
AC_OUTPUT

########################## SUMMARY ############################

AS_IF([test "x$have_xerces" == "xno" || test "x$have_cfitsio" == "xno" ],
  [
   AC_MSG_WARN(
     [

                   Building without xerces-c and/or cfitsio.

                *** The gyoto executable will not be built! ***

      Gyoto can still be used as a C++ library and through the Yorick plug-in
        (although without XML support and/or without certain object types).
     ])
  ])
AS_IF([test "x$have_udunits" == "xno"],
  [
   AC_MSG_WARN(
     [

                         Building without udunits

                     *** Unit conversions disabled ***

     ])
  ])
AS_IF([test "x$have_boost" == "xno"],
  [
   AC_MSG_WARN(
     [

                         Building without boost

                   *** Advanced integrators disabled ***

        The advanced integrators require the Boost library (>=1.53)
            and a C++ compiler that supports the c++11 standard.

     ])
  ])
AS_IF([test "x$have_2f1" == "xno"],
  [
   AC_MSG_WARN(
     [

                         Building without 2F1

                   *** Hypergeometric functions unavailable ***

                      Astrobj::Jet requires ARBLIB or AEAE.

     ])
  ])
