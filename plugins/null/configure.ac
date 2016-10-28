# Whenever this file is edited, reprocess it with autoreconf.

AC_PREREQ([2.69])

# Choose a good name: replace `null' with your own name.
# Also change the e-mail address to your own.
AC_INIT([gyoto-null], [0.0.0], [gyoto@sympa.obspm.fr])
# Below, give the name of some file that will stay in you source tree.
AC_CONFIG_SRCDIR([include/GyotoNull.h])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE
LT_INIT

# Checks for programs.
AC_PROG_INSTALL
PKG_PROG_PKG_CONFIG([])
AS_IF([test "x$PKG_CONFIG" == "x"],
  [AC_MSG_ERROR([Please install pkg-config])])


# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

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
AC_SUBST([DYLIB_VAR])
# Used in yorick/stdplug/Makefile
AC_SUBST([SYS])
AC_SUBST([plugin_sfx])

# Read Gyoto's pkg-config file
PKG_CHECK_MODULES([Gyoto],
   [gyoto],
   [AC_SUBST([YORICK],        [`$PKG_CONFIG --variable=YORICK         gyoto`])
    AC_SUBST([PYTHON],        [`$PKG_CONFIG --variable=PYTHON         gyoto`])
    AC_SUBST([SWIG],          [`$PKG_CONFIG --variable=SWIG           gyoto`])
    AC_SUBST([HOME_LORENE],   [`$PKG_CONFIG --variable=HOME_LORENE    gyoto`])
    AC_SUBST([LORENECPPFLAGS],[`$PKG_CONFIG --variable=LORENECPPFLAGS gyoto`])
    AC_SUBST([LORENELDFLAGS], [`$PKG_CONFIG --variable=LORENELDFLAGS  gyoto`])
    AC_SUBST([soverdir],      [`$PKG_CONFIG --variable=GYOTO_PLUGDIR  gyoto`])
    AC_SUBST([GYOTO],         [`$PKG_CONFIG --variable=GYOTO          gyoto`])
    AC_SUBST([GYOTO_PREFIX],  [`$PKG_CONFIG --variable=prefix         gyoto`])
    AC_SUBST([GYOTO_EPREFIX], [`$PKG_CONFIG --variable=exec_prefix    gyoto`])
    AC_SUBST([Y_INST_HOME],   [`$PKG_CONFIG --variable=Y_INST_HOME    gyoto`])
    AC_SUBST([Y_INST_SITE],   [`$PKG_CONFIG --variable=Y_INST_SITE    gyoto`])
    gyoto_MPICXX=`$PKG_CONFIG --variable=MPICXX gyoto`
    gyoto_CXX=`$PKG_CONFIG --variable=CXX gyoto`
   ],
   [AC_MSG_ERROR([gyoto.pc not found or not functional! 
                  Check PKG_CONFIG_PATH or call the guru])]
)

AC_PROG_CXX(["$gyoto_MPICXX" "$gyoto_CXX" g++ c++ gpp aCC CC cxx cc++ cl.exe FCC KCC RCC xlC_r xlC])

mySUBDIRS="lib"
AC_CONFIG_FILES([Makefile lib/Makefile bin/gyoto],
  [test -e bin/gyoto && chmod a+x bin/gyoto])
AC_SUBST([mySUBDIRS])
AC_OUTPUT
