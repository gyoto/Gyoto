SHELL=/bin/sh
DIR=lib bin
CLEAN_DIR=$(addprefix clean-,$(DIR))
INSTALL_DIR=$(addprefix install-,$(DIR))
UNINSTALL_DIR=$(addprefix uninstall-,$(DIR))

include local_settings

all: $(DIR)

$(DIR) doc:
	cd $@; $(MAKE)

yorick:
	cd $@; yorick -batch make.i; $(MAKE)
check-yorick:
	cd yorick; \
	$(DYLIB_VAR)=../lib:$$$(DYLIB_VAR) \
	GYOTO_CHECK_NODISPLAY=true \
	yorick -i check.i

install: $(INSTALL_DIR)
$(INSTALL_DIR) install-yorick:
	cd $(subst install-,,$@) ; $(MAKE) install

uninstall: $(UNINSTALL_DIR)
$(UNINSTALL_DIR):
	cd $(subst uninstall-,,$@) ; $(MAKE) uninstall

clean: $(CLEAN_DIR) clean-yorick clean-doc
	rm -f example-*.fits
$(CLEAN_DIR) clean-doc:
	cd $(subst clean-,,$@); $(MAKE) clean
clean-yorick:
	cd yorick; yorick -batch make.i; make clean

CHECK_CMD:=$(DYLIB_VAR)=lib:$$$(DYLIB_VAR) bin/gyoto
check:
	$(CHECK_CMD) doc/examples/example-fixed-star.xml \
	   \!example-fixed-star.fits
	$(CHECK_CMD) doc/examples/example-moving-star.xml \
	   \!example-moving-star.fits
	$(CHECK_CMD) doc/examples/example-thin-infinite-disk-BL.xml \
	   \!example-thin-infinite-disk-BL.fits
	$(CHECK_CMD) doc/examples/example-thin-infinite-disk-KS.xml \
	   \!example-thin-infinite-disk-KS.fits
	$(CHECK_CMD) doc/examples/example-torus.xml \
	   \!example-torus.fits

.PHONY: $(DIR) clean $(CLEAN_DIR) $(INSTALL_DIR) yorick clean-yorick install-yorick doc clean-doc install all

