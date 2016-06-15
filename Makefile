# $Id: Makefile,v 1.2 2002-02-28 17:17:43 cbilderback Exp $
# Makefile for the UGM program
# Last modified by $Author: cbilderback $ on $Date: 2002-02-28 17:17:43 $

prefix       = /home/cbilder
host_os      = linux
srcdir       = .
top_srcdir   = .
enable_debug = no

# Set up the include paths
INCPATHS = -I$(prefix)/include -I$(prefix)/include/tiff -I$(prefix)/include/geotiff 
LIBDIRS  = -L$(prefix)/lib

# Libraries we need to link in
LIBS =  -lProjectionMesh -lMathLib -lProjectionIO -lImageLib  -lgeotiff -ltiff -lProjection  -lgctpc -lMiscUtils -lACE -lminipvfs

#SlaveLibs
SLIBS = -lProjectionMesh  -lMiscUtils -lMathLib -lProjectionIO -lImageLib  -lgeotiff -ltiff -lProjection  -lgctpc -lminipvfs

# Linker flags
LDFLAGS   = $(LIBDIRS)
LOADLIBES = $(LIBS)

# Set our compiler options
ifeq ($(enable_debug),yes)
DEBUG = -g -Wall -DTNT_NO_BOUNDS_CHECK 
else
DEBUG = -O3 -Wall -DTNT_NO_BOUNDS_CHECK
#-march=pentiumpro -mcpu=pentiumpro -fomit-frame-pointer -mieee-fp -fschedule-insns2 -finline-functions -frerun-loop-opt -fstrength-reduce -ffast-math -funroll-loops -fexpensive-optimizations -fthread-jumps
endif

# Compiler and other defs
CC   = mpicc
CXX  = mpiCC
CXXFLAGS = $(DEBUG) $(INCPATHS)

# Suffix rules
.SUFFIXES: .o .cpp
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

# Dependencies for the master program
OBJS = Projector.o ProjectionParams.o mastermain.o ProjectorException.o \
       MpiProjector.o BaseProgress.o CLineProgress.o ProjUtil.o Stitcher.o \
       StitcherNode.o inparms.o PVFSProjector.o

SOBJ = Projector.o ProjectionParams.o slavemain.o ProjectorException.o \
       MpiProjectorSlave.o BaseProgress.o ProjUtil.o

all: master slave

master : $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o master $(LIBDIRS) $(LIBS)
slave : $(SOBJ)
	$(CXX) $(CXXFLAGS) $(SOBJ) -o slave $(LIBDIRS) $(SLIBS)


clean:
	rm -f $(OBJS) $(SOBJ) *~ master slave











