# these values filled in by    yorick -batch make.i
Y_MAKEDIR=
Y_EXE=
Y_EXE_PKGS=
Y_EXE_HOME=
Y_EXE_SITE=
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)

# ------------------------------------------------ macros for this package

include ../local_settings

YGYOTO_C = $(shell ls *.C)
GYOTO_C = $(shell ls ../lib/*.C)
GYOTO_OBJ = $(GYOTO_C:.C=.o)
GYOTO_H = $(shell ls ../include/*.h)

PKG_NAME=gyoto
PKG_I=gyoto.i \
	gyoto_Photon.i gyoto_Scenery.i
#gyoto_Kerr.i gyoto_Star.i gyoto_FixedStar.i gyoto_PolishDoughnut.i

OBJS=gyoto_utils.o gyoto_Metric.o gyoto_Astrobj.o gyoto_Photon.o \
	gyoto_Scenery.o gyoto_Screen.o gyoto_Idx.o gyoto_Spectrum.o \
	gyoto_Spectrometer.o gyoto_ThinDisk.o

# gyoto_KerrBL.o gyoto_KerrKS.o gyoto_Star.o gyoto_FixedStar.o gyoto_PolishDoughnut.o gyoto_ThinInfiniteDiskBL.o gyoto_ThinInfiniteDiskKS.o 

# change to give the executable a name other than yorick
PKG_EXENAME=yorick

# PKG_DEPLIBS=-Lsomedir -lsomelib   for dependencies of this package
PKG_DEPLIBS=-L../lib -lgyoto $(LDLIBS)
# set compiler (or rarely loader) flags specific to this package
PKG_CFLAGS= -I../include
PKG_LDFLAGS= -L../lib

# list of additional package names you want in PKG_EXENAME
# (typically Y_EXE_PKGS should be first here)
EXTRA_PKGS=$(Y_EXE_PKGS)

# list of additional files for clean
PKG_CLEAN=

# autoload file for this package, if any
PKG_I_START=
# non-pkg.i include files for this package, if any
PKG_I_EXTRA=gyoto_constants.i gyotoy.i

# -------------------------------- standard targets and rules (in Makepkg)

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plug-in library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# ------------begin C++ source hacks
# must use C++ to load yorick with this C++ package
# this assumes make default CXX macro points to C++ compiler
CXXFLAGS=$(CFLAGS) $(GYOTO_FLAGS)
LD_DLL=$(CXX) $(LDFLAGS) $(PLUG_SHARED)
LD_EXE=$(CXX) $(LDFLAGS) $(PLUG_EXPORT)

# C++ has no standard file extension, supply default make rule(s)
.C.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<
# --------------end C++ source hacks

# simple example:
#myfunc.o: myapi.h
# more complex example (also consider using PKG_CFLAGS above):
#myfunc.o: myapi.h myfunc.c
#	$(CC) $(CPPFLAGS) $(CFLAGS) -DMY_SWITCH -o $@ -c myfunc.c

all: gyoto_std.so

gyoto_std.so: stdplug/gyoto_std.so
	cp stdplug/gyoto_std.so ./

stdplug/gyoto_std.so: stdplug/*.C
	cd stdplug; yorick -batch make.i; make

../lib/$(LIBGYOTO_FILE): ../lib/*.C
	cd ../lib ; make $(LIBGYOTO_FILE)

../lib/libgyoto.a:
	cd ../lib ; make libgyoto.a

gyoto.so: gyoto.i $(OBJS) ../lib/$(LIBGYOTO_FILE)

install::
	cd stdplug; make install
	$(YNSTALL) ygyoto.h ygyoto_supplier.C $(DEST_Y_SITE)/include
	$(YNSTALL) gyotoy.glade $(DEST_Y_SITE)/glade
	$(YNSTALL) gyotoy.py $(DEST_Y_SITE)/python
	$(YNSTALL) gyotoy $(DEST_Y_SITE)/bin
	-ln $(DEST_Y_SITE)/bin/gyotoy $(DESTDIR)$(PREFIX)/bin/
	$(YNSTALL) gyotoy.1 $(DESTDIR)$(PREFIX)/share/man/man1

uninstall::
	cd stdplug; make uninstall
	rm -f $(DEST_Y_SITE)/include/ygyoto*
	rm -f $(DEST_Y_SITE)/glade/gyotoy.glade
	rm -f $(DEST_Y_SITE)/python/gyotoy.py
	rm -f $(DEST_Y_SITE)/bin/gyotoy
	rm -f $(DESTDIR)$(PREFIX)/bin/gyotoy
	rm -f $(DESTDIR)$(PREFIX)/share/man/man1/gyotoy.1

clean::
	cd stdplug; yorick -batch make.i; make clean
	rm -f gyoto_std.so example-patterndisk.fits.gz example-patterndisk.xml

distclean::
	cd stdplug; yorick -batch make.i; make distclean


.PHONY: standalone test-cppprog test-standalone test test-plugin check

# -------------------------------------------------------- end of Makefile
