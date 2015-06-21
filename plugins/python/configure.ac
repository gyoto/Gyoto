#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([gyoto-python], [0.0.0], [gyoto@sympa.obspm.fr])
AC_CONFIG_SRCDIR([include/GyotoPython.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE
LT_INIT

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
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

PKG_CONFIG_PATH+=../../lib
export PKG_CONFIG_PATH

# Read Gyoto's pkg-config file
PKG_CHECK_MODULES([Gyoto],
   [gyoto],
   [
    AC_SUBST([soverdir],      [`$PKG_CONFIG --variable=GYOTO_PLUGDIR  gyoto`])
    AC_SUBST([GYOTO],         [`$PKG_CONFIG --variable=GYOTO          gyoto`])
    AC_SUBST([GYOTO_PREFIX],  [`$PKG_CONFIG --variable=prefix         gyoto`])
    AC_SUBST([GYOTO_EPREFIX], [`$PKG_CONFIG --variable=exec_prefix    gyoto`])
    AC_SUBST([Y_INST_HOME],   [`$PKG_CONFIG --variable=Y_INST_HOME    gyoto`])
    AC_SUBST([Y_INST_SITE],   [`$PKG_CONFIG --variable=Y_INST_SITE    gyoto`])
    gyoto_MPICXX=`$PKG_CONFIG --variable=MPICXX gyoto`
   ],
   [AC_MSG_ERROR([gyoto.pc not found or not functional! 
                  Check PKG_CONFIG_PATH or call the guru])]
)

AS_IF([test "x$gyoto_MPICXX" != "x"], [CXX="$gyoto_MPICXX"])

PC_INIT
PC_PROG_PYTHON
PC_PYTHON_PROG_PYTHON_CONFIG
PC_PYTHON_CHECK_INCLUDES
PC_PYTHON_CHECK_HEADERS
PC_PYTHON_CHECK_CFLAGS
PC_PYTHON_CHECK_LDFLAGS
PC_PYTHON_CHECK_LIBS

AC_SUBST([pluginname], [`basename $PYTHON`])

mySUBDIRS="lib"
AC_CONFIG_FILES([Makefile lib/Makefile bin/gyoto],
  [test -e bin/gyoto && chmod a+x bin/gyoto])
AC_SUBST([mySUBDIRS])
AC_OUTPUT
