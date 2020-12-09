ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I@top_srcdir@/include $(XERCES_CFLAGS) $(UDUNITS_CFLAGS) $(BOOST_CPPFLAGS) \
	-DPACKAGE_STRING="$(PACKAGE_STRING)"
AM_LDFLAGS  = $(XERCES_LIBS) $(PTHREAD_LIBS) $(LDFLAGS_EXEC)
AM_CXXFLAGS = $(PTHREAD_CFLAGS) $(CXXFLAGS_EXEC)
CLEANFILES =

bin_PROGRAMS  =
dist_man_MANS =

if BUILD_GYOTO
bin_PROGRAMS += gyoto
dist_man_MANS += gyoto.1
endif

if HAVE_MPI
AM_LDFLAGS += $(BOOST_MPI_LDFLAGS) $(BOOST_MPI_LIBS) \
	      $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS)
bin_PROGRAMS += gyoto-mpi-worker.@sovers@
dist_man_MANS += gyoto-mpi-worker.@sovers@.1
CLEANFILES += gyoto-mpi-worker.@sovers@.1
endif

gyoto_SOURCES = gyoto.C
gyoto_LDADD   = @top_builddir@/lib/libgyoto@FEATURES@.la
gyoto_CPPFLAGS = $(AM_CPPFLAGS) $(CFITSIO_CFLAGS)
gyoto_LDFLAGS  = $(AM_LDFLAGS) $(CFITSIO_LIBS) -export-dynamic

gyoto_mpi_worker_@sovers@_SOURCES = gyoto-mpi-worker.C
gyoto_mpi_worker_@sovers@_LDADD   = @top_builddir@/lib/libgyoto@FEATURES@.la
gyoto_mpi_worker_@sovers@_CPPFLAGS = $(AM_CPPFLAGS) $(CFITSIO_CFLAGS)
gyoto_mpi_worker_@sovers@_LDFLAGS  = $(AM_LDFLAGS) $(CFITSIO_LIBS) -export-dynamic

gyoto-mpi-worker.@sovers@.1: gyoto-mpi-worker.1
	cp $^ $@

export PATH := .:$(PATH)

EXAMPLE_DIR = $(top_srcdir)/doc/examples/
CHECK_RESULTS_ALL := $(notdir $(patsubst %.xml,%.fits,$(wildcard $(EXAMPLE_DIR)*.xml)))
CHECK_ROTSTAR_RESULTS := $(filter %rotstar3_1.fits, $(CHECK_RESULTS_ALL))
CHECK_WIP_RESULTS := $(filter %-wip.fits, $(CHECK_RESULTS_ALL))
CHECK_RESULTS_SKIP = $(CHECK_ROTSTAR_RESULTS) $(CHECK_WIP_RESULTS)
if !HAVE_2F1
CHECK_RESULTS_SKIP += example-jet.fits example-torusjet.fits
endif
CHECK_RESULTS=$(filter-out $(CHECK_RESULTS_SKIP), $(CHECK_RESULTS_ALL))

CHECK_MPI_RESULTS = $(CHECK_RESULTS:.fits=-mpi.fits)
CHECK_WIP_MPI_RESULTS = $(CHECK_WIP_RESULTS:.fits=-mpi.fits)
CLEANFILES += $(CHECK_RESULTS) $(CHECK_MPI_RESULTS)
CLEANFILES += $(CHECK_WIP_RESULTS) $(CHECK_WIP_MPI_RESULTS)

MPIRUN = mpirun

CHECK_ENV = LD_LIBRARY_PATH=@top_builddir@/lib/.libs/:$(LD_LIBRARY_PATH)
CHECK_RAW_CMD = ./gyoto -pstdplug -r32 -T8 -P0
CHECK_CMD = $(CHECK_ENV) $(CHECK_RAW_CMD)

check-nompi-clean:
	rm -f $(CHECK_RESULTS)

check-mpi-clean:
	rm -f $(CHECK_MPI_RESULTS)

check-wip-nompi-clean:
	rm -f $(CHECK_WIP_RESULTS)

check-wip-mpi-clean:
	rm -f $(CHECK_WIP_MPI_RESULTS)

if BUILD_GYOTO
$(CHECK_RESULTS): gyoto $(EXAMPLE_DIR)$(@:.fits=.xml)
	$(CHECK_CMD) $(EXAMPLE_DIR)$(@:.fits=.xml) \!$@

$(CHECK_MPI_RESULTS): gyoto $(EXAMPLE_DIR)$(@:-mpi.fits=.xml)
	$(CHECK_ENV) $(MPIRUN) $(CHECK_RAW_CMD) -P4 $(EXAMPLE_DIR)$(@:-mpi.fits=.xml) \!$@

$(CHECK_WIP_RESULTS): gyoto $(EXAMPLE_DIR)$(@:.fits=.xml)
	$(CHECK_CMD) $(EXAMPLE_DIR)$(@:.fits=.xml) \!$@

$(CHECK_WIP_MPI_RESULTS): gyoto $(EXAMPLE_DIR)$(@:-mpi.fits=.xml)
	$(CHECK_ENV) $(MPIRUN) $(CHECK_RAW_CMD) -P4 $(EXAMPLE_DIR)$(@:-mpi.fits=.xml) \!$@

check-nompi-compute: $(CHECK_RESULTS)
check-nompi: check-nompi-clean check-nompi-compute

check-mpi-compute: $(CHECK_MPI_RESULTS)
check-mpi: check-mpi-clean check-mpi-compute

check-wip-nompi-compute: $(CHECK_WIP_RESULTS)
check-wip-nompi: check-wip-nompi-clean check-wip-nompi-compute

check-wip-mpi-compute: $(CHECK_WIP_MPI_RESULTS)
check-wip-mpi: check-wip-mpi-clean check-wip-mpi-compute

check: check-nompi
check-wip: check-wip-nompi

if HAVE_LORENE
# make check-lorene check-lorene-mpi
# Will
#   1- copy the Nrotstar directory from HOME_LORENE into .check-lorene;
#   2- configure Nrotstar by copying par_eos.d and par_rot.d from GR/Kepler;
#   3- build Nrotstar/nrotstar;
#   4- run it to build Nrotstar/resu.d;
#   5- ray-trace ../doc/examples/example-*-rotstar3_1.xml twice
#      (with and without MPI).

export HOME_LORENE
.check-lorene/resu.d:
	rm -Rf .check-lorene
	cp -a $(HOME_LORENE)/Codes/Nrotstar ./.check-lorene
	cd .check-lorene && rm -f nrotstar *.o
	cd .check-lorene && cp -f Parameters/GR/Kepler/*.d ./
	cd .check-lorene && $(SED) -i.bak 's/1       graph/0       graph/' par_rot.d
	cd .check-lorene && $(MAKE)
	cd .check-lorene && ./nrotstar

CHECK_ROTSTAR_MPI_RESULTS = $(CHECK_ROTSTAR_RESULTS:.fits=-mpi.fits)
CLEANFILES += $(CHECK_ROTSTAR_RESULTS) $(CHECK_ROTSTAR_MPI_RESULTS)

$(CHECK_ROTSTAR_MPI_RESULTS): gyoto $(EXAMPLE_DIR)$(@:-mpi.fits=.xml) .check-lorene/resu.d
	$(CHECK_ENV) $(MPIRUN) $(CHECK_RAW_CMD) -pstdplug,lorene -P4 $(EXAMPLE_DIR)$(@:-mpi.fits=.xml) \!$@

$(CHECK_ROTSTAR_RESULTS): gyoto $(EXAMPLE_DIR)$(@:.fits=.xml) .check-lorene/resu.d
	$(CHECK_CMD) -pstdplug,lorene -T0 $(EXAMPLE_DIR)$(@:.fits=.xml) \!$@

check-lorene-clean:
	rm -f $(CHECK_ROTSTAR_RESULTS)
check-lorene-compute: $(CHECK_ROTSTAR_RESULTS)
check-lorene: check-lorene-clean check-lorene-compute

check-lorene-mpi-clean:
	rm -f $(CHECK_ROTSTAR_MPI_RESULTS)
check-lorene-mpi-compute: $(CHECK_ROTSTAR_MPI_RESULTS)
check-lorene-mpi: check-lorene-mpi-clean check-lorene-mpi-compute

endif

endif

clean-local:
	-rm -rf .check-lorene

echo-stuff:
	echo $(CHECK_RESULTS_ALL)
	echo $(CHECK_RESULTS)

.PHONY: check-nompi-compute check-nompi-clean check-nompi check-mpi-clean check-mpi-compute check-mpi check check-lorene-clean check-lorene-compute check-lorene check-lorene-mpi-clean check-lorene-mpi-compute check-lorene-mpi check
