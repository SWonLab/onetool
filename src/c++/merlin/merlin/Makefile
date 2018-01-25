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

  TARGET_NAME = Core part of merlin
  TARGET      =
  TARGETS     = libmerlin.a
  TESTS       = runall merlin   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  #HEADERS     = lib_include.h f2c.h kahan.h definition.h \
  #              error.h output_streams.h \
  #              string_util.h string_tokenizer.h \
  #              cephes.h functions.h log_double.h gen_mem.h \
  #              matrix.h fmatrix.h matfunctions.h \
  #              trifmatrix.h sinfo.h histogram.h corinfo.h

  GLOBALS     = AssociationAnalysis.cpp FastAssociation.cpp \
				 AnalysisTask.cpp Conquer.cpp \
				 ConquerHaplotyping.cpp DiseaseModel.cpp \
				 ParametricLikelihood.cpp GenotypeInference.cpp \
				 Houdini.cpp \
				 KongAndCox.cpp Manners.cpp MerlinBitSet.cpp \
				 MerlinCluster.cpp MerlinCore.cpp Merlin.cpp \
				 MerlinError.cpp MerlinFamily.cpp MerlinIBD.cpp \
				 InformationContent.cpp MerlinCache.cpp MerlinModel.cpp \
				 MerlinKinship.cpp MerlinKinship15.cpp \
				 MerlinHaplotype.cpp MerlinMatrix.cpp MerlinParameters.cpp \
				 MerlinPDF.cpp MerlinSimulator.cpp MerlinSimwalk2.cpp \
				 MerlinSort.cpp NPL-ASP.cpp NPL-QTL.cpp \
				 Magic.cpp Mantra.cpp Parametric.cpp QtlModel.cpp \
				 Tree.cpp \
				 TreeBasics.cpp TreeIndex.cpp TreeManager.cpp \
				 TreeInfo.cpp TreeFlips.cpp VarianceComponents.cpp

  SRCS        = $(GLOBALS)

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS     = ${SRCS:.cpp=.h,}

#======================================================================
#   Target: libmerlin.a
#----------------------------------------------------------------------

  libmerlin.a.NAME     = Core part of merlin
  libmerlin.a.TYPE     = LIB
  libmerlin.a.OBJS     = $(OBJS)
  libmerlin.a.CXXFLAGS = -I../../include/merlin/clusters/ -I../../include/merlin/merlin/ -I../../include/merlin/libsrc/ -I../../include/merlin/clusters/ -mssse3
  libmerlin.a.CP       = ../../lib/libmerlin.a

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
