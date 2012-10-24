ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -I@top_builddir@/include $(XERCESCPPFLAGS) $(UDUNITS_CPPFLAGS)
AM_LDFLAGS  = $(XERCESLDFLAGS) $(PTHREAD_LIBS)
AM_CXXFLAGS = -rdynamic $(PTHREAD_CFLAGS)
CLEANFILES=example-*.fits

bin_PROGRAMS  = gyoto
gyoto_SOURCES = gyoto.C
gyoto_LDADD   = @top_builddir@/lib/libgyoto.la
gyoto_CPPFLAGS = $(AM_CPPFLAGS) $(CFITSIOCPPFLAGS)
gyoto_LDFLAGS  = $(AM_LDFLAGS) $(CFITSIOLDFLAGS) -export-dynamic

CHECK_CMD=./gyoto
check: gyoto
	$(CHECK_CMD) ../doc/examples/example-thin-disk.xml \
	   \!example-thin-disk.fits
	$(CHECK_CMD) ../doc/examples/example-complex-astrobj.xml \
	   \!example-complex-astrobj.fits
	$(CHECK_CMD) ../doc/examples/example-fixed-star.xml \
	   \!example-fixed-star.fits
	$(CHECK_CMD) --resolution=32 ../doc/examples/example-moving-star.xml \
	   \!example-moving-star.fits
	$(CHECK_CMD) ../doc/examples/example-page-thorne-disk-BL.xml \
	   \!example-page-thorne-disk-BL.fits
	$(CHECK_CMD) ../doc/examples/example-page-thorne-disk-KS.xml \
	   \!example-page-thorne-disk-KS.fits
	$(CHECK_CMD) ../doc/examples/example-torus.xml \
	   \!example-torus.fits
	$(CHECK_CMD) ../doc/examples/example-polish-doughnut.xml \
	   \!example-polish-doughnut.fits
