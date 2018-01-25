#==========================================================================
# Makefile for oneTool
#
#  Author: Yeunjoo E. Song (yeunjoo.song@gmail.com)
#
#  History: Initial implementation                        - yes Nov 24 2014
#
#  Notes:   Makefile for oneTool application
#
#  Copyright (c) 2014  Sungho Won
#    All Rights Reserved
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

  ONETOOLROOT=../../..

  include $(ONETOOLROOT)/config/Global.make

  TARGET_NAME = "Regression analysis"
  TARGET      =
  TARGETS     = libsegreg.a
  TESTTARGETS = libsegreg.a #test_segreg$(EXE)
  TARPREFIX   = ONETOOL
  TESTS       = #true

#==========================================================================
# Source/object file lists                                                |
#--------------------------------------------------------------------------

  NEW         = types/TypeDescription.cpp

  MODELS      = sub_models/type_sub_model.cpp \
                sub_models/cov_sub_model.cpp \
                sub_models/freq_sub_model.cpp \
                sub_models/transmission_sub_model.cpp \
                sub_models/fpmm_sub_model.cpp \
                sub_models/onset_sub_model.cpp \
                sub_models/resid_sub_model.cpp \
                sub_models/ascertainment_sub_model.cpp \
                sub_models/prev_sub_model.cpp \
                model.cpp

  DATA_MGRS   = PedigreeDataSet.cpp

  CALC        = member_calculator.cpp \
                polygenic_transition_calculator.cpp \
                RegPenetranceCommonAdjustments.cpp \
                RegPenetranceContext.cpp \
                RegPenetranceCalculator.cpp \
                SL_calculator.cpp \
                regressive_peeler.cpp \
                FPMM_peeler.cpp \
                mlm_corr_verifier.cpp \
                mlm_resid_corr_calculator.cpp \
                mlm_peeler.cpp \
                mean_split.cpp \
                freq_estimates.cpp \
                pen_function.cpp \
                segreg_calculator.cpp \
                analysis_iterative_models.cpp

  CONTROL     = parser.cpp \
                analysis_out.cpp \
                AnalysisController.cpp

  SRCS        = $(NEW) $(MODELS) $(DATA_MGRS) $(CALC) $(CONTROL)

  OBJS        = ${SRCS:.cpp=.o}

  HEADERS = ${SRCS:.cpp=.h,}

#======================================================================
#   Target: libsegreg.a                                              |
#----------------------------------------------------------------------

   libsegreg.a.NAME     = "SEGREG library"
   libsegreg.a.TYPE     = LIB
   libsegreg.a.OBJS     = $(OBJS)
   libsegreg.a.CP       = ../../lib/libsegreg.a

#======================================================================
#   Target: test_segreg
#----------------------------------------------------------------------

   test_segreg$(EXE).NAME     = Test of segreg calculators
   test_segreg$(EXE).TYPE     = C++
   test_segreg$(EXE).OBJS     = test_segreg.o
   test_segreg$(EXE).DEP      = libsegreg.a
   test_segreg$(EXE).LDLIBS   = -lsegreg  $(LIB_ALL) 

include $(ONETOOLROOT)/config/Rules.make


# DO NOT DELETE THIS LINE -- make depend depends on it.


-include $(SRCS:%.cpp=$(BUILDDIR)/%.d)
-include $(DEP_SRCS:%.cpp=$(BUILDDIR)/%.d)

