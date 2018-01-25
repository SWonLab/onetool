#==========================================================================
# Makefile for oneTool
#
#  Author:  Yeunjoo E. Song (yeunjoo.song@gmail.com)
#
#  History: Initial implementation                        - yes Nov 24 2014
#
#  Notes:   Makefile for oneTool application
#
#  Copyright (c) 2014  Sungho Won
#    All Rights Reserved
#==========================================================================

#==========================================================================
# General target options
#--------------------------------------------------------------------------
#
# -- These options apply to the entire Makefile, not any single target
#    o TARGET  is used for making a specific target
#    o TARGETS is the list of targets to be make by the recursive rule
#    o VERSION and TARPREFIX set parts of the archive name when the
#      tar rule is called.

  ONETOOLROOT=../..
  include $(ONETOOLROOT)/config/Global.make

  TARGET_NAME = Pedigree Information and QC Tool
  TARGET      =
  TARGETS     = libapp.a onetool$(EXE)
  TESTTARGETS = libapp.a onetool$(EXE)
  TARPREFIX   = ONETOOL
  TESTS       = runall onetool   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  HEADERS     = app_base.h   app_version.h   app_wisard.h   app_sage.h   app_merlin.h    onetool.h

  SRCS        = app_base.cpp app_version.cpp app_wisard.cpp app_sage.cpp app_merlin.cpp

  DEP_SRCS    = onetool.cpp
  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: libapp.a
#----------------------------------------------------------------------

  libapp.a.NAME     = Pedigree oneTool Application Library
  libapp.a.TYPE     = LIB
  libapp.a.OBJS     = $(OBJS)
  libapp.a.CP       = ../lib/libapp.a
  libapp.a.CXXFLAGS = -mssse3 -fpermissive
  libapp.a.LDLIBS   = -lsageglobal

#======================================================================
#   Target: onetool                                                      |
#----------------------------------------------------------------------

  onetool$(EXE).PRE      = rm ../lib/liball.a; ar qT ../lib/liball.a ../lib/lib*.a
  onetool$(EXE).NAME     = "Pedigree Information and QC Tool"
  onetool$(EXE).TYPE     = C++
  onetool$(EXE).INSTALL  = yes
  onetool$(EXE).CXXFLAGS = -mssse3 -fpermissive
  onetool$(EXE).OBJS     = onetool.o
  onetool$(EXE).LDLIBS   = -static -lall -lz -Wl,--whole-archive -lpthread -Wl,--no-whole-archive $(LAPACK) -lgfortran $(LDLIBS_R) -lpthread -ldl -zmuldefs
#  onetool$(EXE).LDLIBS   = -lall -Wl,-Bstatic -lz -lgfortran $(LAPACK) -lgfortran -Wl,-Bdynamic $(LDLIBS_R) -lpthread -ldl -zmuldefs

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
