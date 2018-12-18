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

check-mpi check-mpi-clean check-lorene check-lorene-clean check-lorene-mpi check-lorene-mpi-clean:
if BUILD_GYOTO
	cd bin ; $(MAKE) $@
endif
if HAVE_YORICK
	cd yorick ; $(MAKE) $@
endif

.PHONY: $(SUBDIRS) install ChangeLog
