#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([gyoto], [0.0.4], [gyoto@sympa.obspm.fr])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE(gyoto, 0.0.4)
AC_CONFIG_SRCDIR([include/GyotoKerrKS.h])
AC_CONFIG_MACRO_DIR([m4])
AC_LANG(C++)
LT_INIT([dlopen])

VERSINFO=0:0:0
sovers=`echo $VERSINFO | tr : .`
AC_SUBST([VERSINFO])
AC_SUBST([sovers])
AC_SUBST([VERSION], [`cat ${top_srcdir}/VERSION`])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
DOXYGEN=
AC_CHECK_PROGS([DOXYGEN], [doxygen])

# Check dlopenable library suffix
AC_MSG_CHECKING([for suffix for dynamic library files])
AC_SUBST([target_os])
case "x$target_os" in
     xlinux-gnu)
	#plugin_sfx=so
	DYLIB_VAR=LD_LIBRARY_PATH
	SYS=Linux
	;;
     xdarwin*|xmacos)
	#plugin_sfx=dylib
	DYLIB_VAR=DYLD_LIBRARY_PATH
	SYS=Darwin
	;;
     *) AC_MSG_ERROR([unknown target os $target_os]) ;;
esac
plugin_sfx=so
AC_SUBST([plugin_sfx])
AC_MSG_RESULT([$plugin_sfx])
AC_SUBST([SYS])
AC_SUBST([DYLIB_VAR])

# Check for Yorick
yorick_exe=yorick
YORICK=
AC_MSG_CHECKING([whether to use yorick])
AC_ARG_WITH([yorick],
  [AC_HELP_STRING([--with-yorick[[=yorick-executable]]],
     [whether to use yorick and optionnally the local name (possibly with path)
      of the yorick executable. Default: use yorick if it is found in the
      path.])],
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

# Checks for libraries.
AC_CHECK_LIB([dl], [dlopen])

# Checks for header files.
AC_CHECK_HEADERS([float.h stddef.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor getcwd pow sqrt strdup])

pkg_requires=
# Check for Xerces-c
AC_MSG_CHECKING(whether libxerces should be used)
xerces_headers=
xerces_lib=
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
     [xerces_headers=-I${with_xerces_headers}
      AC_MSG_RESULT([${with_xerces_headers}])],
     [AC_MSG_RESULT([no])])
   AC_MSG_CHECKING(whether xerces library path is provided)
   AC_ARG_WITH([xerces-libs],
     [AS_HELP_STRING([--with-xerces-libs=path],
        [location of xercesc library files])],
     [xerces_lib=-I$with_xerces_libs
      AC_MSG_RESULT([$with_xerces_libs])],
     [AC_MSG_RESULT([no])])
   AC_MSG_CHECKING(for libxerces)
   TMPCPPFLAGS=$CPPFLAGS
   CPPFLAGS="$TMPCPPFLAGS $xerces_headers"
   TMPLDFLAGS=$LDFLAGS
   LDFLAGS="$TMPLDFLAGS $xerces_lib"
   TMPLIBS=$LIBS
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
  ],
  []
)
AS_IF([test "x$have_xerces" == "xno"],
      [AS_IF([test "x$with_xerces" == "xyes"],
             [AC_MSG_ERROR([xerces requested but not found])])
       AC_SUBST(use_xerces,"")
      ],
      [AC_SUBST(use_xerces,[-DGYOTO_USE_XERCES])
       XERCESLDFLAGS="$xerces_lib -lxerces-c"
       XERCESCPPFLAGS="$xerces_headers"
       pkg_requires="${pkg_requires} xerces-c"
      ]
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
     [cfitsio_headers=-I${with_cfitsio_headers}
      AC_MSG_RESULT([${with_cfitsio_headers}])],
     [AC_MSG_RESULT([no])])
   AC_MSG_CHECKING(whether cfitsio library path is provided)
   AC_ARG_WITH([cfitsio-libs],
     [AS_HELP_STRING([--with-cfitsio-libs=path],
        [location of cfitsio library files])],
     [cfitsio_lib=-I$with_cfitsio_libs
      AC_MSG_RESULT([$with_cfitsio_libs])],
     [AC_MSG_RESULT([no])])
   AC_MSG_CHECKING(for libcfitsio)
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
  ],
  []
)
AS_IF([test "x$have_cfitsio" == "xno"],
      [AS_IF([test "x$with_cfitsio" == "xyes"],
             [AC_MSG_ERROR([cfitsio requested but not found])])
       AC_SUBST(use_cfitsio,"")
      ],
      [AC_SUBST(use_cfitsio,[-DGYOTO_USE_CFITSIO])
       CFITSIOLDFLAGS="$cfitsio_lib -lcfitsio"
       CFITSIOCPPFLAGS="$cfitsio_headers"
       pkg_requires="${pkg_requires} cfitsio"
      ]
)
AM_CONDITIONAL([HAVE_CFITSIO], [test "x$have_cfitsio" == "xyes"])
AC_SUBST([CFITSIOLDFLAGS])
AC_SUBST([CFITSIOCPPFLAGS])

# Check for LORENE
AC_MSG_CHECKING([for LORENE])
# first check whether LORENE is required and get HOME_LORENE
require_lorene=no
AC_ARG_WITH([lorene],
  [AS_HELP_STRING([--with-lorene[[=/path/to/lorene]]],
     [whether to use the optional LORENE library, and possibly its path.
      By default, LORENE is used if the HOME_LORENE environment variable
      is set.])],
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
  [AC_SUBST([LIB_LAPACK_DEF],
     [`grep LIB_LAPACK $HOME_LORENE/local_settings`])
   AC_SUBST([LIB_GSL_DEF],
     [`grep LIB_GSL $HOME_LORENE/local_settings`])
   AC_SUBST([LIB_PGPLOT_DEF],
     [`grep LIB_PGPLOT $HOME_LORENE/local_settings`])
   AC_SUBST([LIB_LORENE_DEF],
     [`grep LIB_LORENE $HOME_LORENE/local_settings`])
  ]
)


# Finalization: select output files, subdirs etc.
AC_SUBST([pkg_requires])
mySUBDIRS="lib"
AS_IF([test "x$use_xerces" != "x" && test "x$use_cfitsio" != "x"],
  [mySUBDIRS="$mySUBDIRS bin"])
AC_CONFIG_FILES([Makefile
                 bin/Makefile
                 doc/Makefile
                 lib/Makefile
                 lib/automakefile
 		 lib/Makefile.lorene
                 lib/gyoto.pc
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
  ]
)
AS_IF([test "x$DOXYGEN" != "x"], [mySUBDIRS="$mySUBDIRS doc"])
AC_SUBST([mySUBDIRS])
AC_OUTPUT
