#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Gyoto], [0.1.0], [gyoto@sympa.obspm.fr])
gyoto_test_CXXFLAGS=${CXXFLAGS+set}
gyoto_test_CPPFLAGS=${CPPFLAGS+set}
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([include/GyotoKerrKS.h])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
LT_INIT([dlopen disable-static])
AC_CONFIG_HEADERS([config.h include/GyotoConfig.h]) 

VERSINFO=1:0:0
sovers=`echo $VERSINFO | tr : .`
AC_SUBST([VERSINFO])
AC_SUBST([sovers])

# Check for boost
have_boost=yes
BOOST_REQUIRE([1.53.1], [have_boost=no])

AC_DEFINE_UNQUOTED([GYOTO_SOVERS], ["${sovers}"], [Gyoto ABI version])

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
PKG_PROG_PKG_CONFIG([])
DOXYGEN=
AC_CHECK_PROGS([DOXYGEN], [doxygen])

# Provide additional defaults for user variables

AX_APPEND_COMPILE_FLAGS([[-std=c++11]])

AX_APPEND_COMPILE_FLAGS([-rdynamic], [LDFLAGS])

AC_ARG_ENABLE([native],
  [AC_HELP_STRING([--enable-native],
     [append "-march=native" to CXXFLAGS if supported])],
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
     ])
  ]
)

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
AC_SUBST([DYLIB_VAR])
# Used in yorick/stdplug/Makefile
AC_SUBST([SYS])
AC_SUBST([plugin_sfx])

# Check for Yorick
yorick_exe=yorick
YORICK=
AC_MSG_CHECKING([whether to use yorick])
AC_ARG_WITH([yorick],
  [AC_HELP_STRING([--with-yorick[[=yorick-executable]]],
     [force using yorick and optionnally specify which one [default=yorick,
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

# Only if YORICK available, check for install directories
AS_IF([test "x$YORICK" != "x"],
  [AC_MSG_CHECKING([install path for Yorick platform-dependent files])
   Y_EXE_HOME=`echo write,format=\"%s\",Y_HOME\; quit | ${YORICK}| tail -1`
   Y_EXE_SITE=`echo write,format=\"%s\",Y_SITE\; quit | ${YORICK}| tail -1`
   Y_INST_HOME=
   Y_INST_SITE=
   AC_ARG_WITH([yorick-home],
     [AC_HELP_STRING([--with-yorick-home=DIR],
        [architecture-dependent Yorick files [default: the
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
     [AC_HELP_STRING([--with-yorick-site=DIR],
        [architecture-independent Yorick files
	 [default: the yorick directory if it is under ${prefix}, else
	  ${prefix}/lib/yorick]])
     ],
     [Y_INST_SITE=${with_yorick_site}])
   AC_MSG_RESULT([$Y_INST_SITE])
   AC_SUBST([Y_INST_SITE])

  ]
)

# Checks for libraries.
AC_CHECK_LIB([dl], [dlopen])

# Checks for header files.
# We are not using those tests to far
# AC_CHECK_HEADERS([float.h stddef.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

#BOOST_FIND_HEADER([boost/numeric/odeint/stepper/runge_kutta4_classic.hpp])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
# other functions we may need to check for according to autoscan:
# floor getcwd pow sqrt strdup
AC_CHECK_FUNCS([sincos])

pkg_requires=
pkg_cflags=
pkg_libs=

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

# Check for pthreads
use_pthreads=yes
AC_ARG_ENABLE([pthreads],
  [AS_HELP_STRING([--disable-pthreads],
     [disable POSIX threads parallelization])],
  [AS_IF([test "x$enable_pthreads" = "xno"], [use_pthreads=no])])
AS_IF([test "x$use_pthreads" = "xyes"], [AX_PTHREAD])

# Check for UDUNITS
AC_MSG_CHECKING([whether libudunits2 should be used])
udunits_cppflags=
udunits_ldflags=
UDUNITS_CPPFLAGS=
UDUNITS_LDFLAGS=
AC_ARG_WITH(
  [udunits],
  [AS_HELP_STRING( [--with-udunits[[=DIR]]], 
     [Force using udunits2 installed in a given prefix or force disable it.])],
  [
   AS_IF([test "x$withval" == "xno"],
     [AC_MSG_RESULT(no)],
     [test "x$withval" != "xyes"],
     [udunits_cppflags=-I$withval/include
      udunits_ldflags="-L$withval/lib -ludunits2"
      with_udunits=yes
      AC_MSG_RESULT(yes)],
     [AC_MSG_RESULT(yes)])
  ],
  [with_udunits=check
   AC_MSG_RESULT(probe)
  ]
)
have_udunits=no
AS_IF([test "x$with_udunits" != "xno"],
  [
   AC_MSG_CHECKING(whether udunits2 include path is provided)
   AC_ARG_WITH([udunits-inc],
     [AS_HELP_STRING([--with-udunits-inc=DIR],
        [path to base directory with UDUNITS headers])],
     [with_udunits=yes
      udunits_cppflags=-I${withval}
      AC_MSG_RESULT([${withval}])],
     [AC_MSG_RESULT([no])])
   AC_MSG_CHECKING(whether udunits linker flags are provided)
   AC_ARG_WITH([udunits-lib],
     [AS_HELP_STRING([--with-udunits-lib=LDFLAGS],
        [linker flags for udunits])],
     [with_udunits=yes
      udunits_ldflags="${withval}"
      AC_MSG_RESULT([$withval])],
     [AC_MSG_RESULT([no])])

   # Now check whether or not to use pkg-config for udunits
   AS_IF([test "x$udunits_cppflags" == "x" \
       && test "x$udunits_ldflags" == "x" \
       && test "x$PKG_CONFIG"  != "x"],
     [PKG_CHECK_MODULES([UDUNITS],
         [udunits],
	 [pkg_requires="${pkg_requires} udunits"
          have_udunits=yes
	  UDUNITS_CPPFLAGS=${UDUNITS_CFLAGS}
	  UDUNITS_LDFLAGS=${UDUNITS_LIBS}
         ],
	 [AC_MSG_NOTICE([udunits.pc not found])
         ]
      )
     ],
     [])

   # Failing that, check without pkg-config
   AS_IF([test "x$have_udunits" == "xno"],
     [AC_MSG_CHECKING([for udunits (without pkg-config)])
      TMPCPPFLAGS=$CPPFLAGS
      TMPCFLAGS=$CFLAGS
      TMPLDFLAGS=$LDFLAGS
      TMPLIBS=$LIBS
      CPPFLAGS="$TMPCPPFLAGS $udunits_cppflags"
      LDFLAGS="$TMPLDFLAGS $udunits_ldflags"
      LIBS=-ludunits2
      AC_TRY_LINK(
        [#include <udunits2.h>],
        [ut_set_error_message_handler(ut_ignore);],
        [have_udunits=yes
         AC_MSG_RESULT(yes)],
        [AC_MSG_RESULT(no)])
      CPPFLAGS=$TMPCPPFLAGS
      LDFLAGS=$TMPLDFLAGS
      LIBS=$TMPLIBS
      UDUNITS_LDFLAGS="$udunits_ldflags -ludunits2"
      UDUNITS_CPPFLAGS="$udunits_cppflags"
      pkg_cflags="${pkg_cflags} ${UDUNITS_CPPFLAGS}"
      pkg_libs="${pkg_libs} ${UDUNITS_LDFLAGS}"
     ])
  ],
  []
)
AS_IF([test "x$have_udunits" == "xno"],
  [
   AS_IF([test "x$with_udunits" == "xyes"],
     [AC_MSG_ERROR([
 ****************************************
 *** UDUNITS requested but not found! ***
 ****************************************])
     ]
   )
  ],
  [use_udunits=yes
   AC_DEFINE([HAVE_UDUNITS], [1], [Define to 1 if you have libudunits2])]
)
AC_SUBST([UDUNITS_LDFLAGS])
AC_SUBST([UDUNITS_CPPFLAGS])


# Check for Xerces-c
AC_MSG_CHECKING([whether libxerces should be used])
xerces_headers=
xerces_lib=
XERCESCPPFLAGS=
XERCESLDFLAGS=
AC_ARG_WITH(
  [xerces],
  [AS_HELP_STRING( [--with-xerces[[=xerces-prefix]]], 
     [Force using xerces installed in a given prefix or force disable it.])],
  [
   AS_IF([test "x$with_xerces" == "xno"],
     [AC_MSG_RESULT(no)],
     [test "x$with_xerces" != "xyes"],
     [xerces_headers=-I$with_xerces/include
      xerces_lib=-L$with_xerces/lib
      with_xerces=yes
      AC_MSG_RESULT(yes)],
     [AC_MSG_RESULT(yes)])
  ],
  [with_xerces=check
   AC_MSG_RESULT(probe)
  ]
)
have_xerces=no
AS_IF([test "x$with_xerces" != "xno"],
  [
   AC_MSG_CHECKING(whether xerces include path is provided)
   AC_ARG_WITH([xerces-headers],
     [AS_HELP_STRING([--with-xerces-headers=path],
        [location of xercesc include files (without trailing "/xercesc")])],
     [with_xerces=yes
      xerces_headers=-I${with_xerces_headers}
      AC_MSG_RESULT([${with_xerces_headers}])],
     [AC_MSG_RESULT([no])])
   AC_MSG_CHECKING(whether xerces library path is provided)
   AC_ARG_WITH([xerces-libs],
     [AS_HELP_STRING([--with-xerces-libs=path],
        [location of xercesc library files])],
     [with_xerces=yes
      xerces_lib=-L$with_xerces_libs
      AC_MSG_RESULT([$with_xerces_libs])],
     [AC_MSG_RESULT([no])])

   # Now check whether or not to use pkg-config for xerces
   AS_IF([test "x$xerces_headers" == "x" \
       && test "x$xerces_lib" == "x" \
       && test "x$PKG_CONFIG"  != "x"],
     [PKG_CHECK_MODULES([XERCES],
         [xerces-c],
	 [pkg_requires="${pkg_requires} xerces-c"
          have_xerces=yes
	  XERCESCPPFLAGS=${XERCES_CFLAGS}
	  XERCESLDFLAGS=${XERCES_LIBS}
         ],
	 [AC_MSG_NOTICE([xerces-c.pc not found])
         ]
      )
     ],
     [])

   # Failing that, check without pkg-config
   AS_IF([test "x$have_xerces" == "xno"],
     [AC_MSG_CHECKING([for xerces-c (without pkg-config)])
      TMPCPPFLAGS=$CPPFLAGS
      TMPCFLAGS=$CFLAGS
      TMPLDFLAGS=$LDFLAGS
      TMPLIBS=$LIBS
      CPPFLAGS="$TMPCPPFLAGS $xerces_headers"
      LDFLAGS="$TMPLDFLAGS $xerces_lib"
      LIBS=-lxerces-c
      AC_TRY_LINK(
        [#include <xercesc/util/PlatformUtils.hpp>],
        [xercesc::XMLPlatformUtils::Initialize();],
        [have_xerces=yes
         AC_MSG_RESULT(yes)],
        [AC_MSG_RESULT(no)])
      CPPFLAGS=$TMPCPPFLAGS
      LDFLAGS=$TMPLDFLAGS
      LIBS=$TMPLIBS
      XERCESLDFLAGS="$xerces_lib -lxerces-c"
      XERCESCPPFLAGS="$xerces_headers"
      pkg_cflags="${pkg_cflags} ${XERCESCPPFLAGS}"
      pkg_libs="${pkg_libs} ${XERCESLDFLAGS}"
     ])
  ],
  []
)
AS_IF([test "x$have_xerces" == "xno"],
  [
   AS_IF([test "x$with_xerces" == "xyes"],
     [AC_MSG_ERROR([xerces requested but not found])
     ]
   )
  ],
  [use_xerces=yes
   AC_DEFINE([GYOTO_USE_XERCES], [1], [Define to 1 if you have Xerces])]
)
AC_SUBST([XERCESLDFLAGS])
AC_SUBST([XERCESCPPFLAGS])

# Check for CFITSIO
AC_MSG_CHECKING(whether libcfitsio should be used)
cfitsio_headers=
cfitsio_lib=
AC_ARG_WITH(
  [cfitsio],
  [AS_HELP_STRING( [--with-cfitsio[[=cfitsio-prefix]]], 
     [Force using cfitsio installed in a given prefix or force disable it.])],
  [
   AS_IF([test "x$with_cfitsio" == "xno"],
     [AC_MSG_RESULT(no)],
     [test "x$with_cfitsio" != "xyes"],
     [cfitsio_headers=-I$with_cfitsio/include
      cfitsio_lib=-L$with_cfitsio/lib
      with_cfitsio=yes
      AC_MSG_RESULT(yes)],
     [AC_MSG_RESULT(yes)])
  ],
  [with_cfitsio=check
   AC_MSG_RESULT(probe)
  ]
)
have_cfitsio=no
AS_IF([test "x$with_cfitsio" != "xno"],
  [
   AC_MSG_CHECKING(whether cfitsio include path is provided)
   AC_ARG_WITH([cfitsio-headers],
     [AS_HELP_STRING([--with-cfitsio-headers=path],
        [location of cfitsio include files])],
     [with_cfitsio=yes
      cfitsio_headers=-I${with_cfitsio_headers}
      AC_MSG_RESULT([${with_cfitsio_headers}])],
     [AC_MSG_RESULT([no])])
   AC_MSG_CHECKING(whether cfitsio library path is provided)
   AC_ARG_WITH([cfitsio-libs],
     [AS_HELP_STRING([--with-cfitsio-libs=path],
        [location of cfitsio library files])],
     [with_cfitsio=yes
      cfitsio_lib=-L$with_cfitsio_libs
      AC_MSG_RESULT([$with_cfitsio_libs])],
     [AC_MSG_RESULT([no])])

   # Now check whether or not to use pkg-config for cfitsio
   AS_IF([test "x$cfitsio_headers" == "x" \
       && test "x$cfitsio_lib" == "x" \
       && test "x$PKG_CONFIG"  != "x"],
     [PKG_CHECK_MODULES([CFITSIO],
         [cfitsio],
	 [pkg_requires="${pkg_requires} cfitsio"
          have_cfitsio=yes
	  CFITSIOCPPFLAGS=${CFITSIO_CFLAGS}
	  CFITSIOLDFLAGS=${CFITSIO_LIBS}
         ],
	 [AC_MSG_NOTICE([cfitsio.pc not found])
         ]
      )
     ],
     [])

   # Failing that, check without pkg-config
   AS_IF([test "x$have_cfitsio" == "xno"],
     [AC_MSG_CHECKING([for cfitsio (without pkg-config)])
      TMPCPPFLAGS=$CPPFLAGS
      CPPFLAGS="$TMPCPPFLAGS $cfitsio_headers"
      TMPLDFLAGS=$LDFLAGS
      LDFLAGS="$TMPLDFLAGS $cfitsio_lib"
      TMPLIBS=$LIBS
      LIBS=-lcfitsio
      AC_TRY_LINK(
        [#include <fitsio.h>
         #include <string>],
        [fitsfile* fptr=NULL; std::string file=""; int status=0;
         fits_open_file(&fptr, file.c_str(), 0, &status);],
        [have_cfitsio=yes
         AC_MSG_RESULT(yes)],
        [AC_MSG_RESULT(no)])
      CPPFLAGS=$TMPCPPFLAGS
      LDFLAGS=$TMPLDFLAGS
      LIBS=$TMPLIBS
      CFITSIOLDFLAGS="$cfitsio_lib -lcfitsio"
      CFITSIOCPPFLAGS="$cfitsio_headers"
      pkg_cflags="${pkg_cflags} ${CFITSIOCPPFLAGS}"
      pkg_libs="${pkg_libs} ${CFITSIOLDFLAGS}"
     ]
   )
  ],
  []
)
AS_IF([test "x$have_cfitsio" == "xno"],
  [AS_IF([test "x$with_cfitsio" == "xyes"],
     [AC_MSG_ERROR([cfitsio requested but not found])
     ]
   )
  ],
  [use_cfitsio=yes
   AC_DEFINE([GYOTO_USE_CFITSIO],
     [1],
     [Define to 1 if you have cfitsio])
  ]
)
AM_CONDITIONAL([HAVE_CFITSIO], [test "x$have_cfitsio" == "xyes"])
AC_SUBST([CFITSIOLDFLAGS])
AC_SUBST([CFITSIOCPPFLAGS])

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
   LORENELDFLAGS="-L$HOME_LORENE/Lib -llorene -llorenef77"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_LAPACK | tr -d '\n\r'`"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_GSL | tr -d '\n\r'`"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_PGPLOT | tr -d '\n\r'`"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_LORENE | tr -d '\n\r'`"
   LORENELDFLAGS="$LORENELDFLAGS `$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make LIB_CXX | tr -d '\n\r'`"
   LORENECPPFLAGS=`$MAKE_S HOME_LORENE=$HOME_LORENE -f conftest.make INC | tr -d '\n\r'`
   LORENECPPFLAGS="$LORENECPPFLAGS"
   AC_DEFINE([GYOTO_USE_LORENE],
     [1],
     [Define to 1 if you have Lorene])
   rm conftest.make
   AC_SUBST([LORENECPPFLAGS])
   AC_SUBST([LORENELDFLAGS])
   pkg_cflags="${pkg_cflags} ${LORENECPPFLAGS}"
   pkg_libss="${pkg_libss} ${LORENELDFLAGS}"

  ]
)

# disable doc
AC_ARG_ENABLE([doc],
  [AS_HELP_STRING([--disable-doc],
     [disable building documentation])],
  [AS_IF([test "x$enable_doc" == "xyes"],
     [mkdoc=yes],
     [test "x$enable_doc" == "xno"],
     [mkdoc=no],
     [AC_MSG_ERROR([Please set --disable-doc or --enable-doc only to ``yes'' or ``no''])])
  ],
  [mkdoc=yes])

# Finalization: select output files, subdirs etc.
AC_SUBST([pkg_requires])
AC_SUBST([pkg_cflags])
AC_SUBST([pkg_libs])
mySUBDIRS="lib"
AS_IF([test "x$use_xerces" != "x" && test "x$use_cfitsio" != "x"],
  [mySUBDIRS="$mySUBDIRS bin"])
AC_CONFIG_FILES([Makefile
                 bin/Makefile
                 doc/Makefile
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
  ]
)
AS_IF([test "x$DOXYGEN" != "x" && test "x$mkdoc" == "xyes"], [mySUBDIRS="$mySUBDIRS doc"])
AC_SUBST([mySUBDIRS])
AC_OUTPUT
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

     ])
  ])
