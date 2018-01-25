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

  TARGET_NAME = SAGE Genotype Elimination Model
  TARGET      =
  TARGETS     = libgelim.a
  TESTTARGETS = libgelim.a
  TESTS       = true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = valid_parental_genotypes.h \
                inconsistency_handler.h    \
                geno_eliminate.h           \
                ped_imodel_gen.h

  SRCS        = valid_parental_genotypes.cpp \
                inconsistency_handler.cpp    \
                geno_eliminate.cpp           \
                ped_imodel_gen.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: libgelim.a                                                 |
#----------------------------------------------------------------------

   libgelim.a.NAME     = Genotype Elimination Library
   libgelim.a.TYPE     = LIB
   libgelim.a.OBJS     = $(OBJS)
   libgelim.a.CP       = ../../lib/libgelim.a


include $(ONETOOLROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

