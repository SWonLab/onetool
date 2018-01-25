#==========================================================================
# Makefile for infoQC
#
#  Author:  Yeunjoo E. Song (yeunjoo.song@gmail.com)
#
#  History: Initial implementation                        - yes Nov 24 2014
#
#  Notes:
#    The purpose of this file is to add the appropriate pre-compiler
#    define for the current platform.
#
#    It will append to the CXXFLAGS a precompiler define taking the following
#    form:
#
#      Linux:         __LINUX__
#      SunOS/solaris: __SOLARIS__
#      Windows:       __WIN32__
#      Windows:       __MACOS__
#
#  Copyright (c) 2014 Sungho Won
#    All Rights Reserved
#==========================================================================

# Fetch platform:

UNAME := $(shell uname -s)

# Clear PLATFORMDEFINE:

PLATFORMDEFINE = -Dno_platform

# Assign PLATFORMDEFINE:

ifeq ($(UNAME), Linux)

  PLATFORMDEFINE = -D__LINUX__

endif

ifeq ($(UNAME), SunOS)

  PLATFORMDEFINE = -D__SOLARIS__

endif

ifeq ($(UNAME), CYGWIN_NT-5.0)

  PLATFORMDEFINE = -D__WIN32__

endif

ifeq ($(UNAME), Darwin)

  PLATFORMDEFINE = -D__MACOS__

endif

# Append to CXXFLAGS:

CXXFLAGS := ${CXXFLAGS} ${PLATFORMDEFINE}
