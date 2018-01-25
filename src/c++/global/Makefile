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

  TARGET_NAME = Global utility functions for oneTool
  TARGET      =
  TARGETS     = libglobal.a
  TESTTARGETS = libglobal.a
  TARPREFIX   = ONETOOL
  TESTS       = runall global   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  HEADERS     = app_log.h       option.h   worker.h   analysis.h   common.h \
                interactive.h   io.h       Rconn.h 

  SRCS        = app_log.cpp     option.cpp worker.cpp analysis.cpp common.cpp \
                interactive.cpp io.cpp     Rconn.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: libglobal.a
#----------------------------------------------------------------------

  libglobal.a.NAME     = Utility for oneTool
  libglobal.a.TYPE     = LIB
  libglobal.a.OBJS     = $(OBJS)
  libglobal.a.CXXFLAGS = -mssse3 -fpermissive 
  libglobal.a.CP       = ../lib/libglobal.a

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
