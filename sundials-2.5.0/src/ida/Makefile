# -----------------------------------------------------------------
# $Revision: 1.12 $
# $Date: 2009/03/25 23:10:50 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# Makefile for IDA module
#
# src/ida/Makefile.  Generated from Makefile.in by configure.
# -----------------------------------------------------------------

SHELL = /bin/bash



srcdir       = .
builddir     = .
abs_builddir = /home/arun/Desktop/compuGUT-v01_wSUNDIALS/sundials-2.5.0/src/ida
top_builddir = ../../
prefix       = /home/arun/Desktop/compuGUT-v01_wSUNDIALS/v01/SUNDIALS_LIB
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib

INSTALL        = /usr/bin/install -c
INSTALL_LIB    = ${INSTALL}
INSTALL_HEADER = ${INSTALL} -m 644

LIBTOOL      = $(SHELL) $(top_builddir)/libtool
LIBTOOL_DEPS = config/ltmain.sh

CPP      = cc -E
CPPFLAGS = 
CC       = cc
CFLAGS   = -g -O2
LDFLAGS  = 
LIBS     = -lm 

LAPACK_ENABLED = no

top_srcdir = $(srcdir)/../..

INCLUDES = -I$(top_srcdir)/include -I$(top_builddir)/include

LIB_REVISION = 2:0:0

IDA_LIB = libsundials_ida.la

IDA_SRC_FILES = ida.c ida_ic.c ida_io.c ida_direct.c ida_dense.c ida_band.c ida_spils.c ida_spbcgs.c ida_spgmr.c ida_sptfqmr.c ida_bbdpre.c
IDA_BL_SRC_FILES = ida_lapack.c

IDA_OBJ_FILES = $(IDA_SRC_FILES:.c=.o)
IDA_BL_OBJ_FILES = $(IDA_BL_SRC_FILES:.c=.o)

IDA_LIB_FILES = $(IDA_SRC_FILES:.c=.lo)
IDA_BL_LIB_FILES = $(IDA_BL_SRC_FILES:.c=.lo)

SHARED_LIB_FILES = $(top_builddir)/src/sundials/sundials_band.lo \
                   $(top_builddir)/src/sundials/sundials_dense.lo \
                   $(top_builddir)/src/sundials/sundials_direct.lo \
                   $(top_builddir)/src/sundials/sundials_iterative.lo \
                   $(top_builddir)/src/sundials/sundials_spgmr.lo \
                   $(top_builddir)/src/sundials/sundials_spbcgs.lo \
                   $(top_builddir)/src/sundials/sundials_sptfqmr.lo \
                   $(top_builddir)/src/sundials/sundials_math.lo \
                   $(top_builddir)/src/sundials/sundials_nvector.lo

mkinstalldirs = $(SHELL) $(top_srcdir)/config/mkinstalldirs
rminstalldirs = $(SHELL) $(top_srcdir)/config/rminstalldirs

# ----------------------------------------------------------------------------------------------------------------------

all: $(IDA_LIB)

$(IDA_LIB): shared $(IDA_LIB_FILES)
	@if test "X${LAPACK_ENABLED}" = "Xyes"; then \
	  make lib_with_bl;                          \
	else                                         \
	  make lib_without_bl;                       \
	fi

lib_without_bl: shared $(IDA_LIB_FILES)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) -o $(IDA_LIB) $(IDA_LIB_FILES) $(SHARED_LIB_FILES) -rpath $(libdir) $(LDFLAGS) $(LIBS) -version-info $(LIB_REVISION)

lib_with_bl: shared $(IDA_LIB_FILES) $(IDA_BL_LIB_FILES)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) -o $(IDA_LIB) $(IDA_LIB_FILES) $(IDA_BL_LIB_FILES) $(SHARED_LIB_FILES) -rpath $(libdir) $(LDFLAGS) $(LIBS) -version-info $(LIB_REVISION)

install: $(IDA_LIB)
	$(mkinstalldirs) $(includedir)/ida
	$(mkinstalldirs) $(libdir)
	$(LIBTOOL) --mode=install $(INSTALL_LIB) $(IDA_LIB) $(libdir)
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida.h         $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_direct.h  $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_dense.h   $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_band.h    $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_spbcgs.h  $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_spgmr.h   $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_sptfqmr.h $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_bbdpre.h  $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_spils.h   $(includedir)/ida/
	$(INSTALL_HEADER) $(top_srcdir)/src/ida/ida_impl.h        $(includedir)/ida/
	@if test "X${LAPACK_ENABLED}" = "Xyes"; then \
	  $(INSTALL_HEADER) $(top_srcdir)/include/ida/ida_lapack.h  $(includedir)/ida/ ; \
	fi

uninstall:
	$(LIBTOOL) --mode=uninstall rm -f $(libdir)/$(IDA_LIB)
	rm -f $(includedir)/ida/ida.h
	rm -f $(includedir)/ida/ida_direct.h
	rm -f $(includedir)/ida/ida_dense.h
	rm -f $(includedir)/ida/ida_band.h
	rm -f $(includedir)/ida/ida_lapack.h
	rm -f $(includedir)/ida/ida_spbcgs.h
	rm -f $(includedir)/ida/ida_spgmr.h
	rm -f $(includedir)/ida/ida_sptfqmr.h
	rm -f $(includedir)/ida/ida_bbdpre.h
	rm -f $(includedir)/ida/ida_spils.h
	rm -f $(includedir)/ida/ida_impl.h
	$(rminstalldirs) ${includedir}/ida

shared:
	@cd ${top_builddir}/src/sundials ; \
	 ${MAKE} ;                         \
	 cd ${abs_builddir}

clean:
	$(LIBTOOL) --mode=clean rm -f $(IDA_LIB)
	rm -f $(IDA_LIB_FILES)
	rm -f $(IDA_BL_LIB_FILES)
	rm -f $(IDA_OBJ_FILES)
	rm -f $(IDA_BL_OBJ_FILES)

distclean: clean
	rm -f Makefile

ida.lo: $(srcdir)/ida.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida.c
ida_ic.lo: $(srcdir)/ida_ic.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_ic.c
ida_io.lo: $(srcdir)/ida_io.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_io.c
ida_direct.lo: $(srcdir)/ida_direct.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_direct.c
ida_dense.lo: $(srcdir)/ida_dense.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_dense.c
ida_band.lo: $(srcdir)/ida_band.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_band.c
ida_lapack.lo: $(srcdir)/ida_lapack.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_lapack.c
ida_spils.lo: $(srcdir)/ida_spils.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_spils.c
ida_spbcgs.lo: $(srcdir)/ida_spbcgs.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_spbcgs.c
ida_spgmr.lo: $(srcdir)/ida_spgmr.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_spgmr.c
ida_sptfqmr.lo: $(srcdir)/ida_sptfqmr.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_sptfqmr.c
ida_bbdpre.lo: $(srcdir)/ida_bbdpre.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/ida_bbdpre.c

libtool: $(top_builddir)/$(LIBTOOL_DEPS)
	@cd ${top_builddir} ;                 \
	 ${SHELL} ./config.status --recheck ; \
	 cd ${abs_builddir}
