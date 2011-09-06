SHELL=/bin/sh
include ../local_settings

.SUFFIXES : .o .C

SRC = Astrobj.C Factory.C Register.C SmartPointer.C Utils.C Metric.C \
	Worldline.C FocalPlane.C Photon.C Scenery.C WorldlineIntegState.C \
	Error.C Screen.C Spectrum.C Spectrometer.C

STDPLUG_SRC = KerrBL.C KerrKS.C \
	Star.C FixedStar.C Torus.C \
	ThinInfiniteDiskBL.C ThinInfiniteDiskKS.C \
	PowerLawSpectrum.C BlackBodySpectrum.C

LORENEPLUG_SRC= RotStar3_1.C LorenePlug.C

OBJ = $(SRC:.C=.o)
STDPLUG_OBJ = $(STDPLUG_SRC:.C=.o)
LORENEPLUG_OBJ = $(LORENEPLUG_SRC:.C=.o)

PLUGINS = libgyoto-stdplug.$(DYLIB_SFX)

ifdef HOME_LORENE
PLUGINS += libgyoto-lorene.$(DYLIB_SFX)
endif

all: libgyoto.a $(LIBGYOTO_FILE) $(PLUGINS)


list_src:
	echo $(SRC)

.C.o:
	$(CXX)  -c $(GYOTO_FLAGS) $(CXXFLAGS_G) $(GYOTO_INC) $<

libgyoto.a: $(OBJ)
	rm -f $@
	$(STLIB_CMD) $@ $^

$(LIBGYOTO_FILE): $(OBJ)
	g++ -o $@ $(DYLIB_CFLAG)  $(LIB_CXX) $^

libgyoto-stdplug.$(DYLIB_SFX): $(STDPLUG_OBJ) StdPlug.o $(LIBGYOTO_FILE) Metric.o
	g++ -o $@  $(PLUG_FLAGS)  $(LIB_CXX) $(STDPLUG_OBJ) StdPlug.o

$(LORENEPLUG_OBJ): $(LORENEPLUG_SRC)
	make -f Makefile.lorene $@

libgyoto-lorene.$(DYLIB_SFX): $(LORENEPLUG_OBJ) $(LIBGYOTO_FILE)
	make -f Makefile.lorene

doc: doc.C $(SRC)  *.h
	doxygen doxyfile

clean:
	rm -f $(OBJ) $(STDPLUG_OBJ) $(LORENEPLUG_OBJ) StdPlug.o *.$(DYLIB_SFX)*
	rm -fr Doc
	rm -fr *~
	rm -f libgyoto.a
	rm -f ../include/*~
	rm -Rf Gyoto

$(DESTDIR)$(PREFIX)/include $(DESTDIR)$(PREFIX)/lib \
 $(DESTDIR)$(PREFIX)/lib/gyoto/$(SOVERS):
	install -d $@

install: libgyoto.a $(LIBGYOTO_FILE) $(PLUGINS) \
	            $(DESTDIR)$(PREFIX)/include $(DESTDIR)$(PREFIX)/lib \
	            $(DESTDIR)$(PREFIX)/lib/gyoto/$(SOVERS)
	install -m 0755 $(LIBGYOTO_FILE) $(DESTDIR)$(PREFIX)/lib
	install -m 0755 libgyoto-stdplug.$(DYLIB_SFX) $(DESTDIR)$(PREFIX)/lib/gyoto/$(SOVERS)
	-install -m 0755 libgyoto-lorene.$(DYLIB_SFX) $(DESTDIR)$(PREFIX)/lib/gyoto/$(SOVERS)
	install -m 0644 libgyoto.a $(DESTDIR)$(PREFIX)/lib
	install -m 0644 ../include/Gyoto*.h $(DESTDIR)$(PREFIX)/include

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/libgyoto*
	rm -f $(DESTDIR)$(PREFIX)/include/Gyoto*.h

.PHONY: clean install uninstall