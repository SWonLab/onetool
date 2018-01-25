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
  TARGETS     = libutil.a 
  TESTTARGETS = libutil.a
  TESTS       = runall util   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  HEADERS     = cdflib.h     ftp.h      matrix.h  simd.h    test.h  \
                comb.h       genet.h    memmgr.h  timer.h   data.h  \
                integrate.h  permute.h  sort.h    util.h    dcdflib.h \
                log.h        rand.h     stack.h   vector.h  ext.h   \
                logic.h      regex.h    stat.h    vis.h     new_stat.h \
                matrices/blk.h   matrices/idt.h   matrices/sparse.h  matrices/sublk.h \
                matrices/zero.h  matrices/diag.h  matrices/one.h     matrices/std.h  \
                matrices/sym.h

  SRCS        = comb.cpp    dcdflib.cpp  ftp.cpp    integrate.cpp  logic.cpp  \
                memmgr.cpp  rand.cpp     sort.cpp   stat.cpp       util.cpp   \
                vis.cpp     data.cpp     ext.cpp    genet.cpp      log.cpp    \
                matrix.cpp  permute.cpp  regex.cpp  stack.cpp      test.cpp   \
                timer.cpp   vector.cpp   new_stat.cpp \
                matrices/blk.cpp    matrices/diag.cpp    matrices/idt.cpp  \
                matrices/one.cpp    matrices/sparse.cpp  matrices/std.cpp  \
                matrices/sublk.cpp  matrices/sym.cpp     matrices/zero.cpp

  OBJS        = ${SRCS:.cpp=.o}

    #======================================================================
    #   Target: libutil.a
    #----------------------------------------------------------------------

      libutil.a.NAME     = Utility for oneTool
      libutil.a.TYPE     = LIB
      libutil.a.OBJS     = $(OBJS)
      libutil.a.CXXFLAGS = -mssse3 -fpermissive
      libutil.a.CP       = ../lib/libutil.a
      libutil.a.LDLIBS   = -linput

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
