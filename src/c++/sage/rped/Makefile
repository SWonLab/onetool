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

  TARGET_NAME = SAGE General Referenced Pedigrees and related objects
  TARGET      =
  TARGETS     = librped.a
  TESTTARGETS = librped.a test_rped$(EXE)
  TESTS       = runall rped   #true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  HEADERS     = rped.h   loop.h   fped.h  
#                genome_description.h \
#                Invalidator.h        \
#                loop.h               \
#                rped.h               \
#                rpfile.h             \

  SRCS        = rped.cpp rped_loop.cpp fped_obj.cpp
#                genome_description.cpp \
#                loop.cpp               \
#                rped.cpp               \
#                rpfile.cpp             \
#                rpfile_fortran.cpp     \
#                rpfile_delimited.cpp   \

  DEP_SRCS    = test_rped.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: librped.a                                                 |
#----------------------------------------------------------------------

   librped.a.NAME     = New General Referenced Pedigree Library
   librped.a.TYPE     = LIB
   librped.a.OBJS     = $(OBJS)
   librped.a.CP       = ../../lib/librped.a

#======================================================================
#   Target: test_rped                                                 |
#----------------------------------------------------------------------

   test_rped.NAME     = Test of Rped
   test_rped.TYPE     = C++
   test_rped.OBJS     = test_rped.o 
   test_rped.DEP      = librped.a
   test_rped.LDLIBS   = -lrped


include $(ONETOOLROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

