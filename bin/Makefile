SHELL=/bin/sh

include ../local_settings

.SUFFIXES : .o .C

SRC = gyoto.C
OBJ = $(SRC:.C=.o)

LIB_SRC = $(shell ls ../lib/*.C)
HDRS = $(shell ls ../include/*.h)

EXE_ALL = gyoto
# the gyoto executable requires XERCES
ifeq (,$(findstring -DGYOTO_USE_XERCES,$(GYOTO_FLAGS)))
EXE = $(filter-out gyoto,$(EXE_ALL))
else
EXE = $(EXE_ALL)
endif

all: $(EXE)

.C.o:
	$(CXX)  -c $(GYOTO_FLAGS) -I../include $(CPPFLAGS) $(CXXFLAGS) $<

gyoto: gyoto.o ../lib/$(LIBGYOTO_FILE)
	$(CXX) -o $@ $< -L../lib -lgyoto -ldl $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ)
	rm -f $(EXE_ALL) 
	rm -fr Doc
	rm -fr *~
	rm -f libgyoto.a
	rm -fr gyoto.dSYM

$(DESTDIR)$(PREFIX)/share/man/man1 $(DESTDIR)$(PREFIX)/bin:
	install -d $@

ifeq (,$(findstring -DGYOTO_USE_XERCES,$(GYOTO_FLAGS)))
install:
else
install: $(EXE) $(DESTDIR)$(PREFIX)/bin $(DESTDIR)$(PREFIX)/share/man/man1 gyoto.1
	install -m 0755 $(EXE) $(DESTDIR)$(PREFIX)/bin
	install -m 0644 gyoto.1 $(DESTDIR)$(PREFIX)/share/man/man1
endif

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/bin/gyoto
	rm -f $(DESTDIR)$(PREFIX)/share/man/man1/gyoto.1

.PHONY: clean install uninstall