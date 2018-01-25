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

  TARGET_NAME = Input functions for oneTool
  TARGET      =
  TARGETS     = libinput.a
  TESTTARGETS = libinput.a
  TARPREFIX   = ONETOOL
  TESTS       = runall input   #true

#==========================================================================
# Source/object file lists
#--------------------------------------------------------------------------

  HEADERS     = app_data.h  bgzf.h   dosage.h   filter.h  impute.h  mach.h  \
                meta.h      plink.h  sampvar.h  sim.h  stream.h  typed.h  vcf.h

  SRCS        = app_data.cpp  bgzf.cpp   dosage.cpp   impute.cpp  mach.cpp  \
                meta.cpp      plink.cpp  sampvar.cpp  sim.cpp     stream.cpp \
                typed.cpp     vcf.cpp

  OBJS        = ${SRCS:.cpp=.o}

#======================================================================
#   Target: libinput.a
#----------------------------------------------------------------------

  libinput.a.NAME     = Input for oneTool
  libinput.a.TYPE     = LIB
  libinput.a.OBJS     = $(OBJS)
  libinput.a.CP       = ../lib/libinput.a
  libinput.a.CXXFLAGS = -mssse3 -fpermissive
  #libinput.a.LDLIBS   = $(LIB_ALL)

include $(ONETOOLROOT)/config/Rules.make

# DO NOT DELETE THIS LINE -- make depend depends on it.

-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)
