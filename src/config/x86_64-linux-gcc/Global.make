#==========================================================================
# Makefile for oneTool
#
#  Author:  Yeunjoo E. Song (yeunjoo.song@gmail.com)
#
#  History: Initial implementation                        - yes Nov 24 2014
#
#  Notes:   Global Compiler options
#
#  Copyright (c) 2014  Sungho Won
#    All Rights Reserved
#==========================================================================

#==========================================================================
# Global Compiler options
#--------------------------------------------------------------------------

  CC          = /usr/bin/gcc
  CXX         = /usr/bin/g++
  CFLAGS      = -DDEBUG_VERBOSE -I.
  CXXFLAGS    = -DDEBUG_VERBOSE -I.
  LDLIBS      = # -lefence
  VERBOSE     =3
  OPTS        = sub
  MAKE        = make -j4
  MAKEDEPEND  = /usr/bin/X11/makedepend -Y. -I../lib -I/usr/include
  RM          = rm
  RM_RF       = -Rf
  CP          = cp
  MV          = mv
  AR          = ar
  AR_CREATE   = ruc
  AR_CXX_CREATE = ruc
  RANLIB      = true

  EXE=
  OBJ=.o
  LIB=.a
  LAPACK= -llapack -lblas -lf2c

  EXCEPTIONS = --no_exceptions --no_rtti

  COVERAGE.CFLAGS      = -a
  COVERAGE.CXXFLAGS    = -a
  COVERAGE.LDFLAGS     = -a

  DEBUG.CFLAGS      = -g
  DEBUG.CXXFLAGS    = -g

  STANDARD.CFLAGS   = -O
  STANDARD.CXXFLAGS = -O

  RELEASE.CFLAGS    = -fPIE -O3 -fomit-frame-pointer -funroll-loops
  RELEASE.CXXFLAGS  = -fPIE -O3 -fomit-frame-pointer -funroll-loops
  RELEASE.LDFLAGS   = -static-libgcc 

#==========================================================================
# Global Compiler options for GCC(64) on Linux/x86_64
#--------------------------------------------------------------------------

  CXXFLAGS   := $(CXXFLAGS)  -I/usr/include/  \
               -I$(ONETOOLROOT)/c++/include -I../include \
               -std=c++11 -Wall -Wno-narrowing -gdwarf-2 -gstrict-dwarf  -fpermissive -mssse3 \
               -D_STL_NO_CONCEPT_CHECKS      \
               -I/usr/local/include

  LDFLAGS    :=-L. -L$(ONETOOLROOT)/c++/lib -L../lib $(LDFLAGS) -L
  LDLIBS     := $(LDLIBS) -lieee

  LIB_PLATSPEC = -lpthread
