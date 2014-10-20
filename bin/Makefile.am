ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -I@top_builddir@/include $(XERCESCPPFLAGS) $(UDUNITS_CPPFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS  = $(XERCESLDFLAGS) $(PTHREAD_LIBS)
AM_CXXFLAGS = $(PTHREAD_CFLAGS)
CLEANFILES=example-*.fits

bin_PROGRAMS  =
dist_man_MANS =

if BUILD_GYOTO
bin_PROGRAMS += gyoto
dist_man_MANS += gyoto.1
endif

if HAVE_MPI
bin_PROGRAMS += gyoto-mpi-worker.@sovers@
dist_man_MANS += gyoto-mpi-worker.@sovers@.1
CLEANFILES += gyoto-mpi-worker.@sovers@.1
endif

gyoto_SOURCES = gyoto.C
gyoto_LDADD   = @top_builddir@/lib/libgyoto.la
gyoto_CPPFLAGS = $(AM_CPPFLAGS) $(CFITSIOCPPFLAGS)
gyoto_LDFLAGS  = $(AM_LDFLAGS) $(CFITSIOLDFLAGS) -export-dynamic

gyoto_mpi_worker_@sovers@_SOURCES = gyoto-mpi-worker.C
gyoto_mpi_worker_@sovers@_LDADD   = @top_builddir@/lib/libgyoto.la -lboost_mpi -lboost_serialization
gyoto_mpi_worker_@sovers@_CPPFLAGS = $(AM_CPPFLAGS) $(CFITSIOCPPFLAGS)
gyoto_mpi_worker_@sovers@_LDFLAGS  = $(AM_LDFLAGS) $(CFITSIOLDFLAGS) -export-dynamic

gyoto-mpi-worker.@sovers@.1: gyoto-mpi-worker.1
	cp $^ $@

export PATH := .:$(PATH)

if BUILD_GYOTO
CHECK_CMD = unset GYOTO_PLUGINS && ./gyoto --resolution=32 --nthreads=8
check: gyoto
	$(CHECK_CMD) ../doc/examples/example-thin-disk.xml \
	   \!example-thin-disk.fits
	$(CHECK_CMD) ../doc/examples/example-thin-disk-KS.xml \
	   \!example-thin-disk-KS.fits
	$(CHECK_CMD) ../doc/examples/example-thin-disk-minkowski-cartesian.xml \
	   \!example-thin-disk-minkowski-cartesian.fits
	$(CHECK_CMD) ../doc/examples/example-thin-disk-minkowski-spherical.xml \
	   \!example-thin-disk-minkowski-spherical.fits
	$(CHECK_CMD) ../doc/examples/example-complex-astrobj.xml \
	   \!example-complex-astrobj.fits
	$(CHECK_CMD) ../doc/examples/example-fixed-star.xml \
	   \!example-fixed-star.fits
if HAVE_UDUNITS
	$(CHECK_CMD) ../doc/examples/example-fixed-star-KS.xml \
	   \!example-fixed-star-KS.fits
endif
	$(CHECK_CMD) ../doc/examples/example-fixed-star-minkowski-cartesian.xml \
	   \!example-fixed-star-minkowski-cartesian.fits
	$(CHECK_CMD) ../doc/examples/example-fixed-star-minkowski-spherical.xml \
	   \!example-fixed-star-minkowski-spherical.fits
	$(CHECK_CMD) ../doc/examples/example-startrace.xml \
	   \!example-startrace.fits
	$(CHECK_CMD) ../doc/examples/example-moving-star.xml \
	   \!example-moving-star.fits
	$(CHECK_CMD) ../doc/examples/example-page-thorne-disk-BL.xml \
	   \!example-page-thorne-disk-BL.fits
	$(CHECK_CMD) ../doc/examples/example-page-thorne-disk-KS.xml \
	   \!example-page-thorne-disk-KS.fits
	$(CHECK_CMD) ../doc/examples/example-page-thorne-disk-BL-with-basis.xml \
	   \!example-page-thorne-disk-BL-with-basis.fits
	$(CHECK_CMD) ../doc/examples/example-torus.xml \
	   \!example-torus.fits
	$(CHECK_CMD) ../doc/examples/example-torus-KS.xml \
	   \!example-torus-KS.fits
if HAVE_UDUNITS
	$(CHECK_CMD) ../doc/examples/example-polish-doughnut.xml \
	   \!example-polish-doughnut.fits
endif
endif
