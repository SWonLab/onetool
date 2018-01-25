#==========================================================================
# Makefile for oneTool
#
#  Author:  Yeunjoo E. Song (yeunjoo.song@gmail.com)
#
#  History: Initial implementation                        - yes Jan 24 2018
#
#  Notes:   Global Compiler options
#           2018/01/25 - under construction
#
#  Copyright (c) 2014  Sungho Won
#    All Rights Reserved
#==========================================================================

#==========================================================================
# Global Compiler options
#--------------------------------------------------------------------------

  CC          = gcc
  CXX         = g++
  CFLAGS      = -DDEBUG_VERBOSE -I.
  CXXFLAGS    = -DDEBUG_VERBOSE -I.
  LDLIBS      = # -lefence
  VERBOSE     =3
  OPTS        = sub
  MAKE        = make -j4
  MAKEDEPEND  = /usr/bin/X11/makedepend -Y. -I../lib -I/usr/include
  RM          = rm
  RM_RF       = -Rf
  CP          = cp
  MV          = mv
  AR          = ar
  AR_CREATE   = ruc
  AR_CXX_CREATE = ruc
  RANLIB      = true

  EXE=.exe
  OBJ=.o
  LIB=.a
  LAPACK= -llapack -lblas -lg2c

  EXCEPTIONS = --no_exceptions --no_rtti

  COVERAGE.CFLAGS      = -a -Wall
  COVERAGE.CXXFLAGS    = -a -Wall
  COVERAGE.LDFLAGS     = -a

  DEBUG.CFLAGS      = -g -Wall
  DEBUG.CXXFLAGS    = -g -Wall

  STANDARD.CFLAGS   = -O2 -Wall
  STANDARD.CXXFLAGS = -O2 -Wall

  RELEASE.CFLAGS    = -O2 -Wall -fomit-frame-pointer
  RELEASE.CXXFLAGS  = -O2 -Wall -fomit-frame-pointer
  #RELEASE.LDFLAGS   = -static-libgcc 
  RELEASE.LDFLAGS   =

#==========================================================================
# Global Compiler options for GCC(64) on CYGWIN_NT/x86_64
#--------------------------------------------------------------------------
