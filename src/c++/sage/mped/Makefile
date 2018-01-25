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

  TARGET_NAME = SAGE MultiPedigree
  TARGET      =
  TARGETS     = libmped.a
  TESTTARGETS = libmped.a test_mped$(EXE)
  TESTS       = runall mped   #true


#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = mp.h mpiterator.h  spbaseiter.h  spiterator.h   \
                mpbase.h sp.h sptypes.h mpfwd.h                 \
                spbase.h spiter.h mp_utilities.h

  SRCS        = spbase.cpp           \
                member_base.cpp      \
                family_base.cpp      \
                subpedigree_base.cpp \
                pedigree_base.cpp    \
                pedigree_builder.cpp \
                spbaseiter.cpp       \
                mpbase.cpp           \
                mp_utilities.cpp

  DEP_SRCS    = test_mped.cpp

  OBJS        = $(SRCS:.cpp=.o)

#======================================================================
#   Target: libmped.a                                                 |
#----------------------------------------------------------------------

   libmped.a.NAME     = SAGE MultiPedigree Library
   libmped.a.TYPE     = LIB
   libmped.a.OBJS     = $(OBJS)
   libmped.a.CP       = ../../lib/libmped.a
   libmped.a.CXXFLAGS = -I.

#======================================================================
#   Target: test_mped                                                 |
#----------------------------------------------------------------------

   test_mped$(EXE).NAME     = SAGE MultiPedigree Test
   test_mped$(EXE).TYPE     = C++
   test_mped$(EXE).OBJS     = test_mped.o
   test_mped$(EXE).LDLIBS   = ../../lib/libmped.a
   test_mped$(EXE).DEP      = libmped.a
   test_mped$(EXE).CXXFLAGS = -L.


include $(ONETOOLROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

