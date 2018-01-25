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

  TARGET_NAME = Libraries for merlin
  TARGET      =
  TARGETS     = libmerlinlib.a
  TESTS       = runall merlinlib   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  #HEADERS     = lib_include.h f2c.h kahan.h definition.h \
  #              error.h output_streams.h \
  #              string_util.h string_tokenizer.h \
  #              cephes.h functions.h log_double.h gen_mem.h \
  #              matrix.h fmatrix.h matfunctions.h \
  #              trifmatrix.h sinfo.h histogram.h corinfo.h

  GLOBALS     = BasicHash.cpp Error.cpp FortranFormat.cpp \
				 GenotypeLists.cpp InputFile.cpp IntArray.cpp Hash.cpp \
				 LongArray.cpp Kinship.cpp KinshipX.cpp MapFunction.cpp \
				 MathCholesky.cpp MathDeriv.cpp MathFloatVector.cpp \
				 MathGenMin.cpp MathGold.cpp MathMatrix.cpp MathStats.cpp \
				 MathNormal.cpp MathSVD.cpp MathVector.cpp \
				 MemoryInfo.cpp MiniDeflate.cpp \
				 Parameters.cpp Pedigree.cpp PedigreeAlleleFreq.cpp PedigreeLoader.cpp \
				 PedigreeDescription.cpp PedigreeFamily.cpp PedigreeGlobals.cpp PedigreeTrim.cpp \
				 PedigreePerson.cpp PedigreeTwin.cpp QuickIndex.cpp Random.cpp Sort.cpp \
				 StringArray.cpp StringBasics.cpp StringMap.cpp \
				 StringHash.cpp TraitTransformations.cpp

  SRCS        = $(GLOBALS)

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS     = ${SRCS:.cpp=.h,}

#======================================================================
#   Target: libmerlinlib.a
#----------------------------------------------------------------------

  libmerlinlib.a.NAME     = Library of merlin
  libmerlinlib.a.TYPE     = LIB
  libmerlinlib.a.OBJS     = $(OBJS)
  libmerlinlib.a.CXXFLAGS = -I../../include/merlin/libsrc/ -mssse3 
  libmerlinlib.a.CP       = ../../lib/libmerlinlib.a

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
