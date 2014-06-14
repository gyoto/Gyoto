ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -I@top_srcdir@/include $(XERCESCPPFLAGS) $(UDUNITS_CPPFLAGS)
AM_LDFLAGS  = $(XERCESLDFLAGS) -export-dynamic $(PTHREAD_LIBS) \
	      $(UDUNITS_LDFLAGS)
AM_CXXFLAGS = $(PTHREAD_CFLAGS) -DGYOTO_PREFIX=\"${prefix}\"

# HEADERS: where they are, where to install them
library_includedir=$(includedir)/Gyoto
library_include_HEADERS=$(top_srcdir)/include/*.h

# MAIN LIBRARY
lib_LTLIBRARIES = libgyoto.la
libgyoto_la_SOURCES =  Astrobj.C Factory.C Register.C SmartPointer.C Utils.C Metric.C	\
	WIP.C Worldline.C Photon.C Scenery.C			\
	WorldlineIntegState.C Error.C Screen.C Spectrum.C		\
	Spectrometer.C ComplexSpectrometer.C UniformSpectrometer.C \
	ThinDisk.C Converters.C Functors.C Hooks.C
libgyoto_la_LIBS = $(XERCESLDFLAGS)
libgyoto_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(VERSINFO)

# STANDARD PLUGIN
soverdir=$(pkglibdir)/@sovers@
sover_LTLIBRARIES = libgyoto-stdplug.la
libgyoto_stdplug_la_SOURCES =  KerrBL.C KerrKS.C Minkowski.C \
	Star.C StarTrace.C FixedStar.C Torus.C \
	PowerLawSpectrum.C BlackBodySpectrum.C \
	ComplexAstrobj.C UniformSphere.C \
	StandardAstrobj.C PageThorneDisk.C \
	ThinDiskPL.C PolishDoughnut.C ThinDiskIronLine.C\
	StdPlug.C
# those files are used only if cfitsio is available
cfitsio_stdplug_sources = PatternDisk.C PatternDiskBB.C DynamicalDisk.C \
	 Disk3D.C DynamicalDisk3D.C DirectionalDisk.C
EXTRA_libgyoto_stdplug_la_SOURCES = $(cfitsio_stdplug_sources)
if HAVE_CFITSIO
libgyoto_stdplug_la_SOURCES  +=	$(cfitsio_stdplug_sources)
endif
libgyoto_stdplug_la_LDFLAGS = -module -export-dynamic $(CFITSIOLDFLAGS) \
				$(AM_LDFLAGS)

# LORENE PLUGIN
if HAVE_LORENE
sover_LTLIBRARIES += libgyoto-lorene.la
endif
libgyoto_lorene_la_SOURCES = RotStar3_1.C NumericalMetricLorene.C LorenePlug.C
libgyoto_lorene_la_LDFLAGS = -module -export-dynamic $(LORENELDFLAGS) $(AM_LDFLAGS)
libgyoto_lorene_la_CPPFLAGS = $(AM_CPPFLAGS) $(LORENECPPFLAGS)

# pkg-config file
pkgconfigdir=$(libdir)/pkgconfig
pkgconfig_DATA=gyoto.pc
