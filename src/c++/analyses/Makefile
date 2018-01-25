#==========================================================================
# Makefile for oneTool
#
#  Author:  Sungyoung Lee (me@lsy.io)
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

  TARGET_NAME = Analyses functions for oneTool
  TARGET      =
  TARGETS     = libanalyses.a
  TESTTARGETS = libanalyses.a
  TARPREFIX   = ONETOOL
  TESTS       = runall analyses  #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  HEADERS     = annot.h       design.h     gemma.h    ld.h        mdr.h    pca.h \
                powercalc.h   qtest.h      score.h    setmgr.h    cnv.h    emai.h \
                fqls.h        grm.h        lod.h      mendel.h    pddt.h   ppp.h \
                raretdt.h     scoreGen.h   tridge.h   corr.h      fam.h    fst.h \
                gstest.h      marker.h     mfqls.h    pharaoh.h   qls.h    regr.h \
                scoreRao.h    window.h     ldprune.h

  SRCS        = annot.cpp     design.cpp   gemma.cpp   ld.cpp      mdr.cpp  pca.cpp \
                powercalc.cpp qtest.cpp    score.cpp   setmgr.cpp  cnv.cpp  emai.cpp \
                fqls.cpp      grm.cpp      lod.cpp     mendel.cpp  pddt.cpp ppp.cpp \
                raretdt.cpp   scoreGen.cpp tridge.cpp  corr.cpp    fam.cpp  fst.cpp \
                gstest.cpp    marker.cpp   mfqls.cpp   pharaoh.cpp qls.cpp  regr.cpp \
                scoreRao.cpp  window.cpp   ldprune.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: libanalyses.a
#----------------------------------------------------------------------

  libanalyses.a.NAME     = Analyses for oneTool
  libanalyses.a.TYPE     = LIB
  libanalyses.a.OBJS     = $(OBJS)
  libanalyses.a.CXXFLAGS = -mssse3 -fpermissive
  libanalyses.a.CP       = ../lib/libanalyses.a
  #libanalyses.a.LDLIBS   = $(LIB_ALL)

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
