ACLOCAL_AMFLAGS = -I m4 
AM_CPPFLAGS = -I@top_srcdir@/include $(XERCESCPPFLAGS) $(LORENECPPFLAGS)
AM_LDFLAGS  = $(Gyoto_LDFLAGS) -module -export-dynamic -avoid-version $(LORENELDFLAGS)
AM_CXXFLAGS = $(Gyoto_CFLAGS) -DGYOTO_PREFIX=\"@prefix@\"

soverdir=@libdir@
sover_LTLIBRARIES = libgyoto-null.la
libgyoto_null_la_SOURCES = Null.C NullPlug.C
