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

all: libgyoto.a libgyoto.$(DYLIB_SFX) $(PLUGINS)


list_src:
	echo $(SRC)

.C.o:
	$(CXX)  -c $(GYOTO_FLAGS) $(CXXFLAGS_G) $(GYOTO_INC) $<

libgyoto.a: $(OBJ)
	rm -f $@
	$(STLIB_CMD) $@ $^

libgyoto.$(DYLIB_SFX): $(OBJ)
	g++ -o $@ $(DYLIB_CFLAG)  $(LIB_CXX) $^

libgyoto-stdplug.$(DYLIB_SFX): $(STDPLUG_OBJ) StdPlug.o libgyoto.$(DYLIB_SFX)
	g++ -o $@  $(PLUG_FLAGS)  $(LIB_CXX) $(STDPLUG_OBJ) StdPlug.o

$(LORENEPLUG_OBJ): $(LORENEPLUG_SRC)
	make -f Makefile.lorene $@

libgyoto-lorene.$(DYLIB_SFX): $(LORENEPLUG_OBJ) libgyoto.$(DYLIB_SFX)
	make -f Makefile.lorene

doc: doc.C $(SRC)  *.h
	doxygen doxyfile

clean:
	rm -f $(OBJ) $(STDPLUG_OBJ) $(LORENEPLUG_OBJ) StdPlug.o *.$(DYLIB_SFX)
	rm -fr Doc
	rm -fr *~
	rm -f libgyoto.a
	rm -f ../include/*~
	rm -Rf Gyoto

$(PREFIX)/include $(PREFIX)/lib:
	install -d $@

install: libgyoto.a libgyoto.$(DYLIB_SFX) $(PLUGINS) \
	            $(PREFIX)/include $(PREFIX)/lib
	install -m 0755 libgyoto.$(DYLIB_SFX) $(PREFIX)/lib
	install -m 0755 libgyoto-stdplug.$(DYLIB_SFX) $(PREFIX)/lib
	-install -m 0755 libgyoto-lorene.$(DYLIB_SFX) $(PREFIX)/lib
	install -m 0644 libgyoto.a $(PREFIX)/lib
	install -m 0644 ../include/Gyoto*.h $(PREFIX)/include

uninstall:
	rm -f $(PREFIX)/lib/libgyoto*
	rm -f $(PREFIX)/include/Gyoto*.h

.PHONY: clean install uninstall