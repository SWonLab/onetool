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

  TARGET_NAME = Regression module for merlin
  TARGET      =
  TARGETS     = libregress.a
  TESTS       = runall regress   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  #HEADERS     = lib_include.h f2c.h kahan.h definition.h \
  #              error.h output_streams.h \
  #              string_util.h string_tokenizer.h \
  #              cephes.h functions.h log_double.h gen_mem.h \
  #              matrix.h fmatrix.h matfunctions.h \
  #              trifmatrix.h sinfo.h histogram.h corinfo.h

  GLOBALS     = AutoFit.cpp Regress.cpp \
				 FancyRegression.cpp RegressAnalysis.cpp \
				 RegressKinship.cpp RegressParameters.cpp

  SRCS        = $(GLOBALS)

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS     = ${SRCS:.cpp=.h,}

#======================================================================
#   Target: libregress.a
#----------------------------------------------------------------------

  libregress.a.NAME     = Regression module of merlin
  libregress.a.TYPE     = LIB
  libregress.a.OBJS     = $(OBJS)
  libregress.a.CXXFLAGS = -I../../include/merlin/regress -I../../include/merlin/libsrc/ -I../../include/merlin/merlin/ -mssse3
  libregress.a.CP       = ../../lib/libregress.a

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
