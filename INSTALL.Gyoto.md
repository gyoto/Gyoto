INSTALLING GYOTO
================

0- DISTRIBUTION PACKAGES
========================

Gyoto comes prepackaged for some systems. Although the development
version may be much more advanced at a given time, using pre-compiled
binaries is the recommended way of using Gyoto for most users.

Debian GNU/Linux:

Gyoto is part of Debian since Wheezy (Debian 7.0). Updated packages
are made available through the official backports infrastructure:
  https://backports.debian.org/
Occasionally, even more advanced packages may be available at:
  https://people.debian.org/~thibaut/

On Debian GNU/Linux (and derivatives), try:
 sudo apt-get install gyoto gyoto-doc yorick-gyoto
Developers may be interested in gyoto-dbg and libgyoto0-dev in addition.

On an Apple Mac, you can install MacPorts (http://www.macports.org/) and run:
 sudo port install Gyoto

If Gyoto is not packaged for your system or if you prefer to build
from source, read on.

If the source code was fetched from git, you may
need to first update de timestamp of a few files by running the script
git-post-merge:
   $ ./git-post-merge
See documentation therein.

Then the build process is in a nutshell:
./configure
make
make check
sudo make install

I- INSTALL THE DEPENDENCIES
===========================
   - a C++ compiler. GCC 4.9 works very well. Several features require
     the C++11 standard. Clang/LLVM is discouraged because it does not
     support fenv.h, which sometimes leads to spurious SIFPE in the
     Yorick plug-in:
       https://llvm.org/bugs/show_bug.cgi?id=23707
   - xercesc-3 (recommended, required for the executable):
       http://xerces.apache.org/xerces-c/
   - cfitsio   (required for the executable and some of the Astrobj):
       http://heasarc.gsfc.nasa.gov/fitsio/
   - libudunits2 (recommended, support for unit conversions):
       http://www.unidata.ucar.edu/software/udunits/
   - boost >= 1.53 (required, contains the integrators)
       http://www.boost.org/
   - an MPI implementation (tested with openmpi, optional). MPI uses
     boost features from boost.mpi, you must use the same version as
     boost.mpi is linked to.
   - Yorick    (optional, provides an interface to the Yorick
     interpreted language, allowing to write Gyoto scripts):
       http://yorick.sourceforge.net/
     Yorick users will also need the yorick-yutils add-on
     (https://github.com/frigaut/yorick-yutils) and may need to install
     the yorick-dev package (in particulat Debian/Ubuntu users).
   - Python (optional, provides an interface to the Python interpreted
     language, allowing to write Gyoto scripts). Python 2.7 and 3.4
     have been tested. For building the Python bindings, the Python
     development files are naturally required (sometimes found in the
     python-dev or python3-dev package), as well as NumPy and Swig-2.0:
       https://www.python.org/
       http://www.numpy.org/
       http://www.swig.org/
     Note that although fairly complete, the Python interface is
     likely to change in future releases. Be ready to adapt your
     scripts, or contact us is stability of the API is important for
     you.
   - LORENE (optional, the libgyoto-lorene plug-in can be built later)
       http://www.lorene.obspm.fr/
     On some systems, LORENE must be built with -fPIC (GYOTO as well,
     but this is the default).
   - developers may need the GNU autotools: autoconf, automake, libtool.

II- CONFIGURE GYOTO
===================
If all the dependencies are installed in standard places (/usr or
/usr/local) and if the default prefix (/usr/local) is OK for you, this
should do:
   $ ./configure

You may need to pass some options or configuration variables. To list
the available options:
 ./configure --help
The standard GNU INSTALL file is appended to this file and documents
the most standard and obscure features.

The --enable-release option is important to package maintainers:
without it, the library name will contain "-unreleased". This is to
allow users to compile new versions without overriding the
system-provided library by default. Binaries distributed e.g. by
package managers should be compiled with --enable-release. Also, the
configure script will append flags to the SONAME when features are not
available. This limits the probability of linking with the wrong
version of the library at run time.

To select a different compiler than the default on your system, set
the CC and CXX environment variables accordingly during the configure
step:
   $ CC=gcc-4.8 CXX=g++-4.8 ./configure

Example: assume you want to install in ${HOME}/mysoft, that LORENE is
in ${HOME}/mysoft/Lorene (but HOME_LORENE is not set), and Xerces and
CFITIO are in /opt/local:
   $ ./configure --prefix=${HOME}/mysoft \
                 --with-lorene=${HOME}/mysoft/Lorene \
                 CPPFLAGS=-I/opt/local/include \
                 LDFLAGS=-L/opt/local/lib

If compiling a release (rather than the master branch from git) for
binary distribution, you should use the --enable-release argument to
get rid of "-unreleased" in the SONAME. Don't do that when compiling
anything which is not an official release, though. The unqualified
name libgyoto.* is reserved for full-featured, official releases. The
configure script takes care of adding suffixes when some features are
disabled, e.g. libgyoto-nompi.*. Note that there is no guarantee that
two -unreleased builds are ABI compatible, even if they share the same
SONAME, because the version information is incremented only
immediately before making and official release.

III- BUILD
==========
   $ make
   $ make -C python

IV- TEST
========
Several example files are provided in doc/examples. All the files
without "rotstar3_1" should work out of the box. The rotstar3_1
require the lorene plug-in and some work to create lorene data
files. You can ray-trace all these sceneries (may take up to a couple
of minutes on recent hardware) with:

   $ make check

Don't worry too much for the "severe" warnings.

You can now open the resulting FITS files with ds9
(http://hea-www.harvard.edu/RD/ds9/) or spydr
(http://www.maumae.net/yorick/doc/spydr_intro.php) or any other
FITS-aware image viewer:

   $ spydr example-*.fits

In addition, if Gyoto has been configured to build the Yorick plug-in,
this plug-in is also tested.

V- INSTALL
==========
If installing to a system location (i.e. if you don't
have right access to PREFIX), you need to gain root privileges using,
for instance, su or sudo:
1- using su:
   $ su -
    (type root password)
   # make install
   # make -C python install

2- using sudo:
   $ sudo make install
    (type your password)

Under Linux, if installing to a system location, you may need to also
run
   # ldconfig -v
or
   $ sudo ldconfig -v

VII- Setting your environment
=============================
If installing in a non-standard place (e.g. under your home
directory), you do not need to run ldconfig, but you need to adapt
your environment for instance by adding the following lines to
${HOME}/.profile (replace <gyoto-prefix> by the actual Gyoto
prefix!). One Gyoto file installed in each directory is listed as a
comment:

export PREFIX=<gyoto-prefix>
export PATH=${PREFIX}/bin:${PATH}                       # gyoto
export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH} # libgyoto.so.*
export MANPATH=${PREFIX}/share/man:${MANPATH}           # gyoto.1
export PKG_CONFIG_PATH=${PREFIX}/lib/pkgconfig          # gyoto.pc

Under Mac OS X, LD_LIBRARY_PATH is replaced by DYLD_LIBRARY_PATH.

It goes beyond the scope of this manual to teach you how to set
environment variables; if in doubt ask the local guru or google...

By default, the Yorick plug-in is also installed under ${prefix}. If
Yorick itself is un ${prefix}, then the plug-in will be installed
directly with the rest of Yorick and hence will be found by Yorick. On
the other hand, if Yorick is not under ${prefix}, the plug-in may not
be found immediately by Yorick. Assuming you used the default prefix
(/usr/local), it should be sufficient to create a file named
${HOME}/Yorick/custom.i containing the three following lines:

require, "pathfun.i";
add_y_home,"/usr/local/lib/yorick/";
command_line= process_argv();

Under Debian and Ubuntu GNU/Linux, /usr/local/lib/yorick/ is by
default in Yorick search paths.


VII- More on GNU ./configure script:
====================================

See the generic INSTALL file.
