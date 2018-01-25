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

  TARGET_NAME = "Lodlink"
  TARGET      =
  TARGETS     = liblodlink.a
  TESTTARGETS = liblodlink.a
  TARPREFIX   = ONETOOL
  TESTS       = #true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = mle_sub_model.h    output.h             instructions.h \
                definitions.h      trans_calculator.h   geno_elim.h \
                peeler.h           likelihood.h         max_opt.h \
                results.h          tasks.h \
                genotype_results.h    genotypes.h \
                linkage_results.h     linkage_tests.h \
                homogeneity_results.h homogeneity_tests.h \
                lods.h                analysis.h

  SRCS        = mle_sub_model.cpp output.cpp           instructions.cpp \
                definitions.cpp   trans_calculator.cpp geno_elim.cpp \
                peeler.cpp        likelihood.cpp       max_opt.cpp \
                results.cpp       tasks.cpp \
                genotype_results.cpp    genotypes.cpp \
                linkage_results.cpp     linkage_tests.cpp \
                homogeneity_results.cpp homogeneity_tests.cpp \
                lods.cpp                analysis.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: liblodlink.a
#----------------------------------------------------------------------

   liblodlink.a.NAME     = "LODLINK Library"
   liblodlink.a.TYPE     = LIB
   liblodlink.a.OBJS     = $(OBJS)
   liblodlink.a.CP       = ../../lib/liblodlink.a

#======================================================================
#   Target: test_lodlink
#----------------------------------------------------------------------

   test_lodlink$(EXE).NAME     = "LODLINK test Tool"
   test_lodlink$(EXE).TYPE     = C++
   test_lodlink$(EXE).INSTALL  = yes
   test_lodlink$(EXE).OBJS     = test_lodlink.o
   test_lodlink$(EXE).LDLIBS   = -llodlink $(LIB_SAGE)

include $(ONETOOLROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
