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

  TARGET_NAME = MultiPedigree Correlations with Standard Errors
  TARGET      =
  TARGETS     = libfcor.a
  TESTTARGETS = libfcor.a
  TARPREFIX   = ONETOOL
  TESTS       = #true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = reldefs.h  reltype.h    reltypename.h \
                relcomp.h  relmatrix.h  definitions.h \
                setweight.h \
                pairsetcountinfo.h  pairsetdata.h \
                corcal.h            subcalM.h \
                accal.h \
                output.h            analysis.h

  SRCS        = relstringbank.cpp    reltype.cpp \
                reltypename.cpp      relmatrix.cpp \
                setweight.cpp \
                pairsetcountinfo.cpp pairsetdata.cpp \
                corcal.cpp           subcalM.cpp \
                accal.cpp \
                fcor_output.cpp      fcor_analysis.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: libfcor.a
#----------------------------------------------------------------------

   libfcor.a.NAME     = "FCOR Library"
   libfcor.a.TYPE     = LIB
   libfcor.a.OBJS     = $(OBJS)
   libfcor.a.CP       = ../../lib/libfcor.a

#======================================================================
#   Target: test_fcor
#----------------------------------------------------------------------

   test_fcor$(EXE).NAME     = "FCOR test Tool"
   test_fcor$(EXE).TYPE     = C++
   test_fcor$(EXE).INSTALL  = yes
   test_fcor$(EXE).OBJS     = test_fcor.o
   test_fcor$(EXE).LDLIBS   = -lfcor $(LIB_SAGE)

include $(ONETOOLROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
