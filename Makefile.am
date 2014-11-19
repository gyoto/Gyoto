ACLOCAL_AMFLAGS = -I m4 
SUBDIRS=@mySUBDIRS@
$(SUBDIRS):
	cd $@ ; $(MAKE)

ChangeLog:
	git log --no-merges > ChangeLog

check-mpi check-mpi-clean:
if BUILD_GYOTO
	cd bin ; $(MAKE) $@
endif
if HAVE_YORICK
	cd yorick ; $(MAKE) $@
endif

.PHONY: $(SUBDIRS) install ChangeLog
