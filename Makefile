SHELL=/bin/sh
DIR=lib bin
CLEAN_DIR=$(addprefix clean-,$(DIR))
INSTALL_DIR=$(addprefix install-,$(DIR))
UNINSTALL_DIR=$(addprefix uninstall-,$(DIR))

include local_settings

GYOTO_VERSION=$(subst gyoto-,,`basename $(CURDIR)`)

all: $(DIR) yorick

$(DIR) doc:
	cd $@; $(MAKE)

install: $(INSTALL_DIR) install-yorick
$(INSTALL_DIR):
	cd $(subst install-,,$@) ; $(MAKE) install

uninstall: $(UNINSTALL_DIR) uninstall-yorick
$(UNINSTALL_DIR):
	cd $(subst uninstall-,,$@) ; $(MAKE) uninstall

clean: $(CLEAN_DIR) clean-doc clean-yorick
	rm -f example-*.fits
$(CLEAN_DIR) clean-doc:
	cd $(subst clean-,,$@); $(MAKE) clean

bin/gyoto: bin/*.C
	make bin

lib/$(LIBGYOTO_FILE): lib/*.C
	make lib

CHECK_CMD:=$(DYLIB_VAR)=lib:$$$(DYLIB_VAR) bin/gyoto
check: lib/$(LIBGYOTO_FILE) bin/gyoto check-yorick
	$(CHECK_CMD) doc/examples/example-fixed-star.xml \
	   \!example-fixed-star.fits
	$(CHECK_CMD) --resolution=32 doc/examples/example-moving-star.xml \
	   \!example-moving-star.fits
	$(CHECK_CMD) doc/examples/example-thin-infinite-disk-BL.xml \
	   \!example-thin-infinite-disk-BL.fits
	$(CHECK_CMD) doc/examples/example-thin-infinite-disk-KS.xml \
	   \!example-thin-infinite-disk-KS.fits
	$(CHECK_CMD) doc/examples/example-torus.xml \
	   \!example-torus.fits

yorick yorick/gyoto.so yorick/gyoto_std.so: yorick/*.C yorick/stdplug/*.C
	if which $(YORICK); then cd yorick; $(YORICK) -batch make.i; $(MAKE); fi

check-yorick: yorick/gyoto.so yorick/gyoto_std.so
	if which $(YORICK); then \
	cd yorick; \
	$(DYLIB_VAR)=../lib:$$$(DYLIB_VAR) \
	GYOTO_CHECK_NODISPLAY=true \
	$(YORICK) -i check.i; \
	fi

install-yorick uninstall-yorick: yorick/gyoto.so yorick/gyoto_std.so
	if which $(YORICK); then \
	cd yorick; \
	$(MAKE) $(subst -yorick,,$@); \
	fi

clean-yorick distclean-yorick:
	if which $(YORICK); then \
	cd yorick; \
	$(YORICK) -batch make.i; \
	$(MAKE) $(subst -yorick,,$@); \
	fi

distclean: clean distclean-yorick

tbz: distclean
	cd ..; \
	tar cvf - \
	--exclude=.git --exclude=.gitignore --exclude=.pc --exclude=debian \
	`basename $(CURDIR)` | \
	bzip2 > gyoto_$(GYOTO_VERSION).tbz

.PHONY: all $(DIR) clean $(CLEAN_DIR) install $(INSTALL_DIR) uninstall $(UNINSTALL_DIR) distclean check yorick check-yorick clean-yorick install-yorick uninstall-yorick doc clean-doc distclean-yorick tbz

