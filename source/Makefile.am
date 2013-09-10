## Process this file with automake to produce Makefile.in

EXTRA_DIST = Makefile.am

bin_PROGRAMS = para2

if LOCALLINK
INCLUDES = -I../include -I../../../nucleus -I../../../ajax -I../../../plplot
else
INCLUDES = -I../include -I${embprefix}/include -I${embprefix}/include/eplplot
endif


if ISSHARED
if ISAIXIA64  #???
if LOCALLINK
AM_CFLAGS = -Wl,-bdynamic -Wl,-brtl -L../../../plplot/.libs -L../../../ajax/.libs -L../../../nucleus/.libs -lnucleus -lajax -leplplot -lgmp -lgmpxx
AM_CXXFLAGS = ${AM_CFLAGS}
else
AM_CFLAGS = -Wl,-bdynamic -Wl,-brtl -L${embprefix}/lib -lnucleus -lajax -leplplot -lgmp -lgmpxx
AM_CXXFLAGS = ${AM_CFLAGS}
endif
else
AM_CFLAGS =
AM_CXXFLAGS =
endif
endif


para2_SOURCES = para2.cc alignment.cc


if LOCALLINK
LDADD = ../../../nucleus/libnucleus.la ../../../ajax/libajax.la ../../../plplot/libeplplot.la -lgmp -lgmpxx $(XLIB)
else
LDADD = -L${embprefix}/lib -lnucleus -lajax -leplplot -lgmp -lgmpxx $(XLIB)
endif
