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

  TARGET_NAME = Pedigree Information and QC Tool
  TARGET      =
  TARGETS     = libpedinfo.a
  TESTTARGETS = libpedinfo.a
  TARPREFIX   = ONETOOL
  TESTS       = #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  HEADERS     = relpair.h   loop.h   stats.h   stats_view.h

  SRCS        = relpair.cpp loop.cpp stats.cpp stats_view.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: libpedinfo.a
#----------------------------------------------------------------------

  libpedinfo.a.NAME     = "Pedigree Statistics Library"
  libpedinfo.a.TYPE     = LIB
  libpedinfo.a.OBJS     = $(OBJS)
  libpedinfo.a.CP       = ../../lib/libpedinfo.a

#======================================================================
#   Target: test_pedigree
#----------------------------------------------------------------------

  test_pedigree$(EXE).NAME     = "Pedigree info test Tool"
  test_pedigree$(EXE).TYPE     = C++
  test_pedigree$(EXE).INSTALL  = yes
  test_pedigree$(EXE).OBJS     = test_pedigree.o
  test_pedigree$(EXE).LDLIBS   = -lpedinfo $(LIB_SAGE)


include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
