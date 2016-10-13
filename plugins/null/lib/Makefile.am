# Whenever this file is edited, reprocess it with autoreconf, then
# re-run the configure script.

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I@top_srcdir@/include $(XERCESCPPFLAGS) $(LORENECPPFLAGS)
AM_LDFLAGS  = $(Gyoto_LDFLAGS) -module -export-dynamic -shared -avoid-version $(LORENELDFLAGS)
AM_CXXFLAGS = $(Gyoto_CFLAGS) -DGYOTO_PREFIX=\"@prefix@\"

soverdir=@libdir@

# replace all occurrences of null below by the actual name of the
# plug-in
sover_LTLIBRARIES = libgyoto-null.la
libgyoto_null_la_SOURCES = Null.C InitPlug.C
