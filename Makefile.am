ACLOCAL_AMFLAGS = -I m4 
SUBDIRS=@mySUBDIRS@

@VERBATIM@ifneq ($(wildcard plugins/python/Makefile),)
@VERBATIM@    SUBDIRS += plugins/python
@VERBATIM@endif

list-dirs:
	echo $(SUBDIRS)

$(SUBDIRS):
	cd $@ ; $(MAKE)

ChangeLog:
	git log --no-merges > ChangeLog

check-mpi check-mpi-clean check-lorene check-lorene-clean check-lorene-mpi check-lorene-mpi-clean check-python check-python-clean check-python-mpi check-python-mpi-clean:
if BUILD_GYOTO
	cd bin ; $(MAKE) $@
endif
if HAVE_PYTHON
	cd python ; $(MAKE) $@
endif
if HAVE_YORICK
	cd yorick ; $(MAKE) $@
endif

.PHONY: $(SUBDIRS) install ChangeLog

# make sure SUBDIRS are built in the right order:

bin python yorick: lib/libgyoto@FEATURES@.la
python: doc/xml/index.xml
lib/libgyoto@FEATURES@.la: lib
doc/xml/index.xml: doc

# build necessary things before calling check targets:
GYOTO_PLUGINS_LA = lib/libgyoto-stdplug.la
if HAVE_PYTHON
GYOTO_PYTHON_LA := plugins/python/lib/$(shell sed -n \
	's/^sover_LTLIBRARIES[[:space:]]*=[[:space:]]*//p' \
	plugins/python/lib/Makefile)
GYOTO_PLUGINS_LA += $(GYOTO_PYTHON_LA)
endif
if HAVE_LORENE
GYOTO_PLUGINS_LA += lib/libgyoto-lorene.la
endif

check check-recursive check-mpi check-lorene check-lorene-mpi check-python check-python-mpi: \
	lib/libgyoto@FEATURES@.la $(GYOTO_PLUGINS_LA)
if HAVE_PYTHON
check check-recursive check-mpi check-python check-python-mpi: \
	python/testbed/stamp
endif
if HAVE_MPI
check-mpi check-lorene-mpi check-python-mpi: \
	bin/gyoto-mpi-worker.@sovers@@EXEEXT@
endif
python/testbed/stamp:
	$(MAKE) -C python testbed
SUBDIR_ARTIFACTS=lib/libgyoto@FEATURES@.la \
	$(GYOTO_PLUGINS_LA) \
	bin/gyoto-mpi-worker.@sovers@@EXEEXT@
$(SUBDIR_ARTIFACTS):
	$(MAKE) -C $(dir $@) $(notdir $@)
