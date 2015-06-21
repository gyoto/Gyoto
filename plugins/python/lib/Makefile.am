ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -I@top_srcdir@/include -DGYOTO_NO_DEPRECATED \
	$(PYTHON_INCLUDES)
AM_LDFLAGS  = $(Gyoto_LDFLAGS) -module -export-dynamic \
	$(PYTHON_LDFLAGS) $(PYTHON_LIBS)
AM_CXXFLAGS = $(Gyoto_CFLAGS) -DGYOTO_PREFIX=\"@prefix@\" \
	$(PYTHON_CFLAGS)


soverdir=@soverdir@
sover_LTLIBRARIES = libgyoto-@pluginname@.la
libgyoto_@pluginname@_la_SOURCES = \
	PythonPlug.C Base.C Spectrum.C Metric.C Standard.C ThinDisk.C

if DEFAULT
install-data-hook:
	cd $(DESTDIR)$(soverdir) && \
	  rm -f libgyoto-python.$(GYOTO_PLUGIN_SFX) && \
	  $(LN_S) libgyoto-$(pluginname).$(GYOTO_PLUGIN_SFX) \
	          libgyoto-python.$(GYOTO_PLUGIN_SFX)
uninstall-hook:
	rm -f $(DESTDIR)$(soverdir)/libgyoto-python.$(GYOTO_PLUGIN_SFX)
endif

# HEADERS: where they are, where to install them
library_includedir=$(GYOTO_PREFIX)/include/Gyoto
library_include_HEADERS=$(top_srcdir)/include/*.h
