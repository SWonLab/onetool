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

  TARGET_NAME = SAGE Multi Locus
  TARGET      =
  TARGETS     = libmlocus.a 
  TESTTARGETS = libmlocus.a test_mlocus$(EXE) 
  TESTS       = runall mlocus

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = genotype.h   phenotype.h   phmodel.h   penmodel.h   imodel.h \
                mfile.h

  SRCS        = genotype.cpp phenotype.cpp phmodel.cpp penmodel.cpp imodel.cpp \
                mfile.cpp

  DEP_SRCS    = test_mlocus.cpp test_phmodel.cpp test_penmodel.cpp test_imodel.cpp test_mfile.cpp

  OBJS        = $(SRCS:.cpp=.o)

#======================================================================
#   Target: libmlocus.a blob                                          |
#----------------------------------------------------------------------

   libmlocus.a.NAME     = SAGE Multi Locus Library
   libmlocus.a.TYPE     = LIB
   libmlocus.a.OBJS     = $(OBJS)
   libmlocus.a.CP       = ../../lib/libmlocus.a

#======================================================================
#   Target: test_mlocus                                               |
#----------------------------------------------------------------------

   test_mlocus$(EXE).NAME     = SAGE Multi Locus Test
   test_mlocus$(EXE).TYPE     = C++
   test_mlocus$(EXE).DEP      = libmlocus.a
   test_mlocus$(EXE).OBJS     = test_mlocus.o test_phmodel.o test_penmodel.o test_imodel.o test_mfile.o
   test_mlocus$(EXE).LDLIBS   = -lmlocus
   test_mlocus$(EXE).CXXFLAGS = -L.

include $(ONETOOLROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

