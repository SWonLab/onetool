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

  TARGET_NAME = Cluster module for merlin
  TARGET      =
  TARGETS     = libclusters.a
  TESTS       = runall clusters   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  #HEADERS     = lib_include.h f2c.h kahan.h definition.h \
  #              error.h output_streams.h \
  #              string_util.h string_tokenizer.h \
  #              cephes.h functions.h log_double.h gen_mem.h \
  #              matrix.h fmatrix.h matfunctions.h \
  #              trifmatrix.h sinfo.h histogram.h corinfo.h

  GLOBALS     = HaploFamily.cpp HaploGraph.cpp \
				 HaploSet.cpp HaploTree.cpp Likelihood.cpp \
				 SparseLikelihood.cpp Unknown.cpp

  SRCS        = $(GLOBALS)

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS     = ${SRCS:.cpp=.h,}

#======================================================================
#   Target: libclusters.a
#----------------------------------------------------------------------

  libclusters.a.NAME     = Cluster module of merlin
  libclusters.a.TYPE     = LIB
  libclusters.a.OBJS     = $(OBJS)
  libclusters.a.CXXFLAGS = -I../../include/merlin/clusters -I../../include/merlin/libsrc/ -I../../include/merlin/merlin/ -mssse3
  libclusters.a.CP       = ../../lib/libclusters.a

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
