# INSTALLING GYOTO

## 0- Installing precompiled packages

Gyoto comes prepackaged for some systems. Although the development
version may be much more advanced at a given time, using pre-compiled
binaries is the recommended way of using Gyoto for most users.

### Debian GNU/Linux

Gyoto is part of Debian since Wheezy (Debian 7.0). Updated packages
are made available through the official [backports infrastructure]
(https://backports.debian.org/). Occasionally, even more advanced
packages may be available at https://people.debian.org/~thibaut/.

You can get a list of available packages with

    apt-cache search gyoto

The most recent packages will install about everything in Gyoto with

    sudo apt-get install gyoto

### Ubuntu

Gyoto is also part of Ubuntu at least since Raring (13.04). Updated
versions are sometimes provided on our [personal package archive
(PPA)] (https://launchpad.net/~paumard/+archive/ubuntu/gyoto/). Check
what version is available there compared to what is available for your
version of Ubuntu and follow instructions on that page to add this PPA
to your system.

You can get a list of available packages with

    apt-cache search gyoto

The most recent packages will install about everything in Gyoto with

    sudo apt-get install gyoto

### Mac OS X

Precompiled binaries (or at least automatic compilation) is provided
through [MacPorts] (http://www.macports.org/). With MacPorts installed, run:

    sudo port sync
    sudo port install Gyoto

To get MPI parallelization, you must first install Boost with one of
its MPI variant, then Gyoto with the same variant:

    sudo port install Boost +openmpi
    sudo port install Gyoto +openmpi

## 1- Building from source: installing the dependencies

If Gyoto is not packaged for your system or if you prefer to build
from source, read on.

The first step is to install the dependencies. Please refer to
[BUGS.md](BUGS.md) for known bugs in some versions of the following
dependencies.

Gyoto requires:

   - a C++ compiler. GCC 4.9 and above work very well. Several
     features require the C++11 standard. Clang/LLVM is discouraged,
     see [BUGS.md](BUGS.md).
   - xercesc-3 (recommended, required for the executable):
       http://xerces.apache.org/xerces-c/
   - cfitsio   (required for the executable and some of the Astrobj):
       http://heasarc.gsfc.nasa.gov/fitsio/
   - libudunits2 (recommended, support for unit conversions):
       http://www.unidata.ucar.edu/software/udunits/
   - boost >= 1.53 (required, contains the integrators). See [BUGS.md](BUGS.md).
       http://www.boost.org/
   - an implementation of the Gauss hypergeometric function 2F1
     (optional, required for all Astrobj using KappaDistributionSynchrotronSpectrum), one of:
     + ARBLIB: http://arblib.org (in that case, compile and install
       ARBLIB and its dependencies and use the --with-arblib*
       configure options);
     + AEAE: http://cpc.cs.qub.ac.uk/summaries/AEAE_v1_0.html (in that
       case unpack the AEAE source code somewhere and use the
       --with-aeae configure option).
   - an MPI implementation (tested with openmpi, optional). MPI uses
     boost features from boost.mpi, you must use the same version as
     boost.mpi is linked to.
   - Yorick    (optional, provides an interface to the Yorick
     interpreted language, allowing to write Gyoto scripts):
       http://yorick.sourceforge.net/
     Yorick users will also need the yorick-yutils add-on
     (https://github.com/frigaut/yorick-yutils) and may need to install
     the yorick-dev package (in particulat Debian/Ubuntu users).
   - Python 3 (optional, provides an interface to the Python
     interpreted language, allowing to write Gyoto scripts). Python
     3.7 and 3.8 have been tested. For building the Python bindings,
     the Python development files are naturally required (sometimes
     found in the python-dev or python3-dev package), as well as NumPy
     and Swig-2.0:
       https://www.python.org/
       http://www.numpy.org/
       http://www.swig.org/
     Note that although fairly complete, the Python interface is
     likely to change in future releases. Be ready to adapt your
     scripts, or contact us is stability of the API is important for
     you.
   - LORENE (optional, the libgyoto-lorene plug-in can be built later):
       https://www.lorene.obspm.fr/
     On some systems, LORENE must be built with -fPIC (GYOTO as well,
     but this is the default).
   - developers may need the GNU autotools: autoconf, automake, libtool.
   - Eigen3 (required for the polarisation):
       https://packages.debian.org/fr/sid/libeigen3-dev

For Debian and its derivatives (incl. Ubuntu), you can install all
those dependencies with:

    sudo apt-get install build-essential yorick-dev yorick-yutils \
    libxerces-c-dev libcfitsio-dev libudunits2-dev libboost-dev \
    libboost-mpi-dev libflint-arb-dev libflint-dev mpi-default-dev \
    python3-dev python3-setuptools swig3.0 python3-numpy python3-matplotlib \
    doxygen pkg-config liblorene-dev lorene-codes-src gfortran g++ libeigen3-dev

## 2- Downloading the source code

The source code is available from
[Github](https://github.com/gyoto/Gyoto):

    git clone https://github.com/gyoto/Gyoto.git

(This obviously requires git to be installed on your system, on Debian
and derivtives use `sudo apt-get install git`).

Then the build process is, in a nutshell, after having installed the
dependencies:

    ./git-post-merge
    ./configure
    make
    sudo make install
    sudo ldconfig

The rest of this file details each step.


## 3- Fixing the timestamps

Unfortunately git does not preserve the timestamps of files, which
confuses the the build system. The easiest way to do that is running a
provided script each time you pull from our repository:

    ./git-post-merge

This script contains instructions to automate this step if you plan of
pulling again from github in the future.

Alternatively, you could recreate the autotools-generated files using
`autoreconf`. This requires the development tools autoconf, automake,
libtool, and is really necessary only for developpers who modified the
the build system (configure.ac, */Makefile.am...)

## 4- Configuring Gyoto

If all the dependencies are installed in standard places (/usr or
/usr/local) and if the default prefix (/usr/local) is OK for you, this
should do:

    ./configure

You may need to pass some options or configuration variables. To list
the available options:

    ./configure --help

ARBLIB is known to be installed under various names depending on the
Linux distribution. If using ARBLIB (see "Installing the dependencies"
above), you may need to set the `--with-arblib-ldflags` variable to the
correct name, e.g.

    ./configure --with-arblib-ldflags=-larb

The standard GNU INSTALL file is provided next to this file and documents
the most standard and obscure features.

The `--enable-release` option is reserved for pre-compiled package
maintainers. In short, don't use it, it is for us alone. Without this
option, the library name will contain "-unreleased". This is to allow
users to compile new versions without overriding the system-provided
library by default. Binaries distributed e.g. by package managers
should be compiled with --enable-release, but only when compiling code
in the `stable` branch. Also, the configure script will append flags
to the SONAME when features are not available,
e.g. libgyoto-nompi.*. This limits the probability of linking with the
wrong version of the library at run time. Note that there is no
guarantee that two -unreleased builds are ABI compatible, even if they
share the same SONAME, because the version information is incremented
only immediately before making and official release.

To select a different compiler than the default on your system, set
the CC and CXX environment variables accordingly during the configure
step:

    CC=gcc-4.8 CXX=g++-4.8 ./configure

Example: assume you want to install in `${HOME}/mysoft`, that LORENE is
in `${HOME}/mysoft/Lorene` (but `HOME_LORENE` is not set), and Xerces and
CFITIO are in `/opt/local`:

    ./configure --prefix=${HOME}/mysoft \
                --with-lorene=${HOME}/mysoft/Lorene \
                CPPFLAGS=-I/opt/local/include \
                LDFLAGS=-L/opt/local/lib

On Debian or Ubuntu, with all the dependencies installed as above,
this should do:

    ./configure --with-arblib \
    --with-lorene=/usr/lib/`dpkg-architecture -qDEB_HOST_MULTIARCH`/lorene

## 5- Building Gyoto

    make

## 6- Testing

Gyoto includes a detailed check suite, including atomic tests writen
in Yorick and Python as well as full ray-tracing tests. MPI tests and
LORENE tests are run using separate Makefile targets. To run all the
tests (which assumes that both Gyoto was configures with both MPI and
LORENE):

    make check check-lorene check-mpi check-lorene-mpi

Don't worry too much for the "severe" warnings.

You can now open the resulting FITS files with ds9
(http://hea-www.harvard.edu/RD/ds9/) or spydr
(http://www.maumae.net/yorick/doc/spydr_intro.php) or any other
FITS-aware image viewer:

    spydr example-*.fits

## 7- Installing

If installing to a system location (i.e. if you don't
have right access to PREFIX), you need to gain root privileges using,
for instance, su or sudo:

Using su:

    su - # (type root password)
    make install
    make -C python install

Using sudo:

    sudo make install
    (type your password)

Under Linux, if installing to a system location, you may need to also
run

    ldconfig -v

as root (so most likely `sudo ldconfig -v`).

## 8- Setting your environment

If installing in a non-standard place (e.g. under your home
directory), you do not need to run ldconfig, but you need to adapt
your environment for instance by adding the following lines to
`${HOME}/.profile` (replace `<gyoto-prefix>` by the actual Gyoto
prefix!). One Gyoto file installed in each directory is listed as a
comment:

    export PREFIX=<gyoto-prefix>
    export PATH=${PREFIX}/bin:${PATH}                       # gyoto
    export LD_LIBRARY_PATH=${PREFIX}/lib:${LD_LIBRARY_PATH} # libgyoto*.so.*
    export MANPATH=${PREFIX}/share/man:${MANPATH}           # gyoto.1
    export PKG_CONFIG_PATH=${PREFIX}/lib/pkgconfig          # gyoto.pc

Under Mac OS X, `LD_LIBRARY_PATH` is replaced by `DYLD_LIBRARY_PATH`.

It goes beyond the scope of this manual to teach you how to set
environment variables; if in doubt ask the local guru or google...

By default, the Yorick plug-in is also installed under `${prefix}`. If
Yorick itself is in `${prefix}`, then the plug-in will be installed
directly with the rest of Yorick and hence will be found by Yorick. On
the other hand, if Yorick is not under `${prefix}`, the plug-in may not
be found immediately by Yorick. Assuming you used the default prefix
(`/usr/local`), it should be sufficient to create a file named
`${HOME}/Yorick/custom.i` containing the three following lines:

    require, "pathfun.i";
    add_y_home,"/usr/local/lib/yorick/";
    command_line= process_argv();

Under Debian and Ubuntu GNU/Linux, `/usr/local/lib/yorick/` is by
default in Yorick search paths.
