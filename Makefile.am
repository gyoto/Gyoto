ACLOCAL_AMFLAGS = -I m4 
SUBDIRS=@mySUBDIRS@
$(SUBDIRS):
	cd $@ ; $(MAKE)

ChangeLog:
	git log --no-merges > ChangeLog

.PHONY: $(SUBDIRS) install ChangeLog
