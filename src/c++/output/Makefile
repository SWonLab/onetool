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

  TARGET_NAME = Output functions for oneTool
  TARGET      =
  TARGETS     = liboutput.a
  TESTTARGETS = liboutput.a
  TARPREFIX   = ONETOOL
  TESTS       = runall output   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  HEADERS     = exporter.h result.h wis.h 

  SRCS        = exporter.cpp result.cpp wis.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: liboutput.a
#----------------------------------------------------------------------

  liboutput.a.NAME     = Output for oneTool
  liboutput.a.TYPE     = LIB
  liboutput.a.OBJS     = $(OBJS)
  liboutput.a.CP       = ../lib/liboutput.a
  liboutput.a.CXXFLAGS = -mssse3 -fpermissive
  liboutput.a.LDLIBS   = -linput

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
