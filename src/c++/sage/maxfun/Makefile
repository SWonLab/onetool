#========================================================================== 
# Makefile for ONETOOL
#
#  Author: Yeunjoo E. Song (yeunjoo.song@gmail.com)
#
#  History: Initial implementation                        - yes Aug 24 2015 
#
#  Copyright (c) 2014  Sungho Won
#==========================================================================

#========================================================================== 
# General target options                                                  |
#--------------------------------------------------------------------------
#
# -- These options apply to the entire Makefile, not any single target
#    o TARGET  is used for making a specific target
#    o TARGETS is the list of targets to be make by the recursive rule
#    o VERSION and TARPREFIX set parts of the archive name when the
#      tar rule is called.

  ONETOOLROOT= ../../..
  include $(ONETOOLROOT)/config/Global.make

  TARGET_NAME = "Maximization of function"
  TARGET      =
  TARGETS     = libmaxfun.a
  TESTTARGETS = libmaxfun.a
  TARPREFIX   = MAXFUN
  TESTS       = runall maxfun

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  SRCS        = maxfun.cpp cmaxfun.cpp maxfun_sub_model.cpp \
                CovarianceMatrix.cpp Parameter.cpp \
                ParameterInput.cpp Submodel.cpp \
                TransformationSubmodel.cpp ParameterMgr.cpp \
                DebugCfg.cpp RunCfg.cpp SequenceCfg.cpp \
                maxfun_Results.cpp APIMaxFunction.cpp Maximizer.cpp

  OBJS        = ${SRCS:.cpp=.o}

  #DEP_SRCS    = test_maxfun.cpp test_maxfun2.cpp

  HEADERS     = ${SRCS:.cpp=.h}

#======================================================================
#   Target: libmaxfun.a                                               |
#----------------------------------------------------------------------

   libmaxfun.a.NAME     = "MAXFUN Library"
   libmaxfun.a.TYPE     = LIB
   libmaxfun.a.OBJS     = ${OBJS}
   libmaxfun.a.CP       = ../../lib/libmaxfun.a
   libmaxfun.a.CXXFLAGS = -static

#======================================================================
#   Target: test_maxfun                                               |
#----------------------------------------------------------------------

   test_maxfun.NAME      = "Simple MAXFUN example"
   test_maxfun.INSTALL   = yes
   test_maxfun.TYPE      = C++
   test_maxfun.SRCS      = test_maxfun.cpp
   test_maxfun.OBJS      = test_maxfun.o
   test_maxfun.DEP       = libmaxfun.a
   test_maxfun.LDLIBS    = $(LIB_ALL) 


#======================================================================
#   Target: test_maxfun2                                                     |
#----------------------------------------------------------------------

   test_maxfun2.NAME      = Simple Maxfun example
   test_maxfun2.INSTALL   = yes
   test_maxfun2.TYPE      = C++
   test_maxfun2.SRCS      = test_maxfun2.cpp
   test_maxfun2.OBJS      = test_maxfun2.o
   test_maxfun2.DEP       = libmaxfun.a
   test_maxfun2.LDLIBS    = $(LIB_ALL) 


include $(ONETOOLROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

