dnl Copyright 2018 Thibaut Paumard
dnl
dnl This file is part of Gyoto.
dnl
dnl Gyoto is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl Gyoto is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.


dnl Define a number of --with options for a given library.
dnl Parameters:
dnl    $1: lower case name of the library. The various "with" options
dnl        will use it: --with-$1, --with-$1-headers and --with-$1-libs
dnl    $2: name of the library as known by pkg_config (think $2.pc)
dnl    $3: linker flag to link with this library (e.g. `-lfoo')
dnl    $4: set of include directives to test for this library (for AC_TRY_LINK)
dnl    $5: main body to test for this library (for AC_TRY_LINK)
dnl Outputs:
dnl    Where $1 is `libname':
dnl    have_libname: yes if found, no otherwise
dnl    with_libname: yes if required, no or probe otherwise
dnl    libname_headers: set of -I flags
dnl    libname_lib: set of -L flags
dnl    LIBNAME_CFLAGS: CPP flags for compiling with this library
dnl    LIBNAME_LIBS: linker flags to link with this library
dnl Side effects:
dnl    also updates pkg_cflags and pkg_libs
AC_DEFUN([GYOTO_ARG_LIB],
[
AC_MSG_CHECKING([whether $1 should be used])
gy_headers=
gy_lib=
gy_CFLAGS=
gy_LIBS=
gy_with=$with_$1

AC_ARG_WITH(
  [$1],
  [AS_HELP_STRING( [--with-$1[[=$1-prefix]]], 
     [Force using $1 installed in a given prefix or force disable it.])],
  [
   AS_IF([test "x$withval" == "xno"],
     [AC_MSG_RESULT(no)],
     [test "x$withval" != "xyes"],
     [gy_headers=-I$withval/include
      gy_lib=-L$withval/lib
      gy_with=yes
      AC_MSG_RESULT(yes)],
     [AC_MSG_RESULT(yes)])
  ],
  [gy_with=check
   AC_MSG_RESULT(probe)
  ]
)

gy_have=no
AS_IF([test "x$gy_with" != "xno"],
  [
   AC_MSG_CHECKING(whether $1 include path is provided)
   AC_ARG_WITH([$1-headers],
     [AS_HELP_STRING([--with-$1-headers=path],
        [location of $1 include files as a colon-separated path ('--with-$1-headers=/path1:...:/pathn') or list of preprocessor flags ('--with-$1-headers="-I/path1 ... -I/pathn")'])],
     [gy_with=yes
      gy_headers=${withval/:/ -I}
      case "x$gy_headers" in
      	   x\-I*) ;;
	   x) ;;
	   x*) gy_headers=-I$gy_headers ;;
      esac
      AC_MSG_RESULT([${gy_headers}])],
     [AC_MSG_RESULT([no])])
   AC_MSG_CHECKING(whether $1 library path is provided)
   AC_ARG_WITH([$1-libs],
     [AS_HELP_STRING([--with-$1-libs=path],
        [location of $1 library files as a colon-separated path ('--with-$1-libs=/path1:...:/pathn') or list of preprocessor flags ('--with-$1-libs="-L/path1 ... -L/pathn")'])],
     [gy_with=yes
      gy_lib=${withval/:/ -L}
      case "x$gy_lib" in
      	   x\-L*) ;;
	   x) ;;
	   x*) gy_lib=-L$gy_lib ;;
      esac
      AC_MSG_RESULT([$gy_lib])],
     [AC_MSG_RESULT([no])])

   # Now check whether or not to use pkg-config for this library
   AS_IF([test "x$gy_headers" == "x" \
       && test "x$gy_lib" == "x" \
       && test "x$PKG_CONFIG"  != "x"],
     [PKG_CHECK_MODULES(translit([$1], [a-z], [A-Z]),
         [$2],
	 [pkg_requires="${pkg_requires} $2"
          gy_have=yes
         ],
	 [AC_MSG_NOTICE([$2.pc not found])
         ]
      )
     ],
     [])

   # Failing that, check without pkg-config
   AS_IF([test "x$gy_have" == "xno"],
     [AC_MSG_CHECKING([for ]translit([$1], [a-z], [A-Z])[ (without pkg-config)])
      TMPCPPFLAGS=$CPPFLAGS
      TMPCFLAGS=$CFLAGS
      TMPLDFLAGS=$LDFLAGS
      TMPLIBS=$LIBS
      CPPFLAGS="$TMPCPPFLAGS $gy_headers"
      LDFLAGS="$TMPLDFLAGS $gy_lib"
      LIBS=$3
      AC_TRY_LINK(
        [$4],
	[$5],
        [gy_have=yes
         $1[]_headers=$gy_headers
         $1[]_lib=$gy_lib
         translit([$1], [a-z], [A-Z])[]_LIBS="$gy_lib $3"
         translit([$1], [a-z], [A-Z])[]_CFLAGS="$gy_headers"
         pkg_cflags="${pkg_cflags} ${translit([$1], [a-z], [A-Z])[]_CFLAGS}"
         pkg_libs="${pkg_libs} ${translit([$1], [a-z], [A-Z])[]_LIBS}"
         AC_MSG_RESULT(yes)],
        [AC_MSG_RESULT(no)])
      CPPFLAGS=$TMPCPPFLAGS
      LDFLAGS=$TMPLDFLAGS
      LIBS=$TMPLIBS
     ])
  ],
  []
)

have_$1=$gy_have
with_$1=$withval

AS_IF([test "x$gy_have" == "xno"],
  [AS_IF([test "x$gy_with" == "xyes"],
     [AC_MSG_ERROR([$1 requested but not found])
     ],
     []
   )
   $1[]_headers=""
   $1[]_lib=""
   translit([$1], [a-z], [A-Z])[]_LIBS=""
   translit([$1], [a-z], [A-Z])[]_CFLAGS=""
  ],
  [use_$1=yes
   AC_DEFINE([GYOTO_USE_]translit([$1], [a-z], [A-Z]),
     [1],
     [Define to 1 if you have $1])
  ]
)
AM_CONDITIONAL([HAVE_]translit([$1], [a-z], [A-Z]), [test "x$have_$1" == "xyes"])
AC_SUBST(translit([$1], [a-z], [A-Z])[_LIBS])
AC_SUBST(translit([$1], [a-z], [A-Z])[_CFLAGS])


])
