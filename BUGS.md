# Known bugs affecting Gyoto

This file is to list known bugs in Gyoto or its dependencies that
significantly affect Gyoto but are not easy to fix or workaround in
Gyoto.

## Partial uninstallation of Python extension

After having installed and uninstalled Gyoto from source with "make
uninstall", empty directories may remain in a Python dist-packages
directory (typically /usr/local/lib/python3.??/dist-packages/), which
causes subsequent installation of Python packages from source to fail
(Gyoto and any other setuptools-based package). Make sure to remove
such empty directories.

                        -- Thibaut Paumard, Mon Dec 09 2024.

## Clang/LLVM compiler (all versions), possibly other compilers:

Gyoto relies on proper delivery of SIGFPE when arithmetic errors occur
(division by zero, floating point overflows...). It does so by
enabling the right exceptions using feenableexcept from GNU fenv.h.

However, some compilers (clang/llvm, in particular) do not support
setting the floating-point environment and will happily optimize code
assuming that it is safe to divide by zero. We recommend staying away
from these compilers until this issue is fixed:
  - https://llvm.org/bugs/show_bug.cgi?id=23707
  - https://llvm.org/bugs/show_bug.cgi?id=8100

Alternatively, you may try to reconfigure Gyoto to not use fenv.h
and/or to not use any optimization: make clean ./configure
 --without-fenv CFLAGS="-g -O0" CXXFLAGS="-g -O0" make
and see whether it works better for you.

If you experience such spurious SIGFPE inside the Yorick plug-in, an
ugly workaround is to call gyoto.fedisableexcept() right before the
offending Gyoto call. If at all possible, consider recompiling Gyoto
as explained above.

Note that loading some plug-ins (e.g. the lorene plug-in) may change
or reset the default floating-point environment. This is the case if
the plug-in, or one of the libraries it links with, has been compiled
with optimizations such as -ffast-math or gfortran.

                        -- Thibaut Paumard, Thu Dec 08 2016.
