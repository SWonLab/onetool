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

  ONETOOLROOT=../../..
  include $(ONETOOLROOT)/config/Global.make

  TARGET_NAME = Global utility functions for SAGE
  TARGET      =
  TARGETS     = libsageglobal.a
  TESTTARGETS = libsageglobal.a test_sinfo$(EXE)
  TESTS       = runall global   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  #HEADERS     = lib_include.h f2c.h kahan.h definition.h \
  #              error.h output_streams.h \
  #              string_util.h string_tokenizer.h \
  #              cephes.h functions.h log_double.h gen_mem.h \
  #              matrix.h fmatrix.h matfunctions.h \
  #              trifmatrix.h sinfo.h histogram.h corinfo.h

  GLOBALS     = error.cpp output_streams.cpp \
                string_util.cpp string_tokenizer.cpp

  NUMERICS    = combinatorial.cpp ndtr.cpp ndtri.cpp gamma.cpp \
                igam.cpp igami.cpp chdtr.cpp fmatrix.cpp \
                corinfo.cpp

  LSFS        = LSFnameMgr.cpp LSFattr.cpp LSFtypes.cpp \
                LSF.cpp LSFfile.cpp aparser.cpp

  SRCS        = $(GLOBALS) $(NUMERICS) $(LSFS)

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS     = ${SRCS:.cpp=.h,}

#======================================================================
#   Target: libsageglobal.a
#----------------------------------------------------------------------

  libsageglobal.a.NAME     = Global Utility for SAGE
  libsageglobal.a.TYPE     = LIB
  libsageglobal.a.OBJS     = $(OBJS)
  libsageglobal.a.CP       = ../../lib/libsageglobal.a

#======================================================================
#   Target: test_sinfo
#----------------------------------------------------------------------

  test_sinfo$(EXE).NAME     = SampleInfo test program
  test_sinfo$(EXE).TYPE     = C++
  test_sinfo$(EXE).OBJS     = test_sinfo.o
  test_sinfo$(EXE).DEP      = libsageglobal.a
  test_sinfo$(EXE).LDLIBS   = -lsageglobal

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
