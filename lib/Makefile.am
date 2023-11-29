ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I@top_srcdir@/include $(XERCES_CFLAGS) $(UDUNITS_CFLAGS) $(BOOST_CPPFLAGS) $(CFITSIO_CFLAGS) $(ARBLIB_CFLAGS) $(AEAE_CFLAGS) $(EIGEN_CFLAGS)
AM_LDFLAGS  = $(XERCES_LIBS) -export-dynamic $(PTHREAD_LIBS) \
	      $(UDUNITS_LIBS) $(CFITSIO_LIBS) $(ARBLIB_LIBS)
AM_CXXFLAGS = $(PTHREAD_CFLAGS) -DGYOTO_PKGLIBDIR=\"${pkglibdir}\"
if HAVE_LOCALPKGLIBDIR
AM_CXXFLAGS += -DGYOTO_LOCALPKGLIBDIR=\"${localpkglibdir}\"
endif
if HAVE_MPI
AM_LDFLAGS += $(BOOST_MPI_LDFLAGS) $(BOOST_MPI_LIBS) \
	      $(BOOST_SERIALIZATION_LDFLAGS) $(BOOST_SERIALIZATION_LIBS)
endif


# HEADERS: where they are, where to install them
library_includedir=$(includedir)/Gyoto
library_include_HEADERS=$(top_srcdir)/include/*.h $(top_builddir)/include/*.h

# MAIN LIBRARY
lib_LTLIBRARIES = libgyoto@FEATURES@.la
libgyoto@FEATURES@_la_SOURCES =  \
	Value.C Property.C Object.C \
	Astrobj.C Factory.C Register.C SmartPointer.C Utils.C Metric.C	\
	WIP.C Worldline.C Photon.C Scenery.C			\
	WorldlineIntegState.C Error.C Screen.C Spectrum.C		\
	Spectrometer.C ComplexSpectrometer.C UniformSpectrometer.C \
	StandardAstrobj.C ThinDisk.C Converters.C Functors.C Hooks.C \
	GridData2D.C FitsRW.C
libgyoto@FEATURES@_la_LIBS = $(XERCES_LIBS)
libgyoto@FEATURES@_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(VERSINFO)

# STANDARD PLUGIN
soverdir=$(pkglibdir)/@sovers@
sover_LTLIBRARIES = libgyoto-stdplug.la
libgyoto_stdplug_la_CPPFLAGS = $(AM_CPPFLAGS) -DGYOTO_PLUGIN=stdplug
libgyoto_stdplug_la_SOURCES =  KerrBL.C KerrKS.C Minkowski.C \
	ChernSimons.C RezzollaZhidenko.C Hayward.C SchwarzschildHarmonic.C \
	Star.C StarTrace.C FixedStar.C InflateStar.C \
	Torus.C OscilTorus.C \
	PowerLawSpectrum.C BlackBodySpectrum.C \
	ThermalBremsstrahlungSpectrum.C \
	ThermalSynchrotronSpectrum.C \
	PowerLawSynchrotronSpectrum.C KappaDistributionSynchrotronSpectrum.C \
	ComplexAstrobj.C UniformSphere.C \
	ComplexMetric.C Shift.C PageThorneDisk.C \
	ThinDiskPL.C PolishDoughnut.C ThinDiskIronLine.C\
	DeformedTorus.C EquatorialHotSpot.C Jet.C Blob.C Plasmoid.C FreeStar.C ThickDisk.C\
	ThinDiskProfile.C SphericalAccretion.C StdPlug.C


# those files need cfitsio for some functionality
cfitsio_stdplug_sources = PatternDisk.C PatternDiskBB.C \
	DynamicalDisk.C DynamicalDiskBolometric.C \
	Disk3D.C DynamicalDisk3D.C DirectionalDisk.C XillverReflection.C \
	FlaredDiskSynchrotron.C ThinDiskGridIntensity.C
EXTRA_libgyoto_stdplug_la_SOURCES = $(cfitsio_stdplug_sources)

libgyoto_stdplug_la_SOURCES  +=	$(cfitsio_stdplug_sources)

libgyoto_stdplug_la_LDFLAGS = -module -export-dynamic $(AM_LDFLAGS) -avoid-version

# LORENE PLUGIN
if HAVE_LORENE
sover_LTLIBRARIES += libgyoto-lorene.la
endif
libgyoto_lorene_la_SOURCES = RotStar3_1.C NumericalMetricLorene.C \
	NeutronStar.C NeutronStarAnalyticEmission.C \
	NeutronStarModelAtmosphere.C LorenePlug.C
libgyoto_lorene_la_LDFLAGS = -module -export-dynamic $(LORENELDFLAGS) $(AM_LDFLAGS) -avoid-version
libgyoto_lorene_la_CPPFLAGS = $(AM_CPPFLAGS) $(LORENECPPFLAGS) -DGYOTO_PLUGIN=lorene

# pkg-config file
pkgconfigdir=$(libdir)/pkgconfig
pkgconfig_DATA=gyoto.pc

uninstall-hook:
	-rmdir $(DESTDIR)$(soverdir)
	-rmdir $(DESTDIR)$(pkglibdir)
	-rmdir $(DESTDIR)$(library_includedir)
