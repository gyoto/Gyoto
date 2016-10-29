# Whenever this file is edited, reprocess it with autoreconf, then
# re-run the configure script.

ACLOCAL_AMFLAGS = -I m4

# Add $(LORENECPPFLAGS) and $(LORENELDFLAGS) or other flags as needed.
# Unless your plug-in needs some deprecated method names, specify
# -DGYOTO_NO_DEPRECATED.
AM_CPPFLAGS = -I@top_srcdir@/include $(XERCESCPPFLAGS)
AM_LDFLAGS  = $(Gyoto_LDFLAGS) -module -export-dynamic -shared -avoid-version
AM_CXXFLAGS = $(Gyoto_CFLAGS) -DGYOTO_NO_DEPRECATED

# Replace all occurrences of null below by the actual name of the
# plug-in.
# List all the C++ source files in libgyoto_<coolname>_la_SOURCES.
sover_LTLIBRARIES = libgyoto-null.la
libgyoto_null_la_SOURCES = Null.C InitPlug.C
