#==========================================================================
# Makefile for merlin parts in oneTool
#
#  Author:  Yeunjoo E. Song (yeunjoo.song@gmail.com)
#
#  History: Initial implementation                        - yes Nov 24 2014
#           gelim added                                   - syl Oct 18 2016
#
#  Notes:   oneTool sage-wide make
#
#  Copyright (c) 2014  Sungho Won
#    All Rights Reserved  
#==========================================================================

VERSION    = 1
PATCHLEVEL = 0
SUBLEVEL   = 0
VERBOSE    = 3
BUILDOPTS  =

ONETOOLROOT=../..
include $(ONETOOLROOT)/config/Global.make

CORE       = merlin libsrc clusters 

TARGETS = $(CORE)

#
# Global rules for INTERNAL build nodes.

#
all: subs

subs:
	@unset MAKEFLAGS; set -e; DIR=`pwd` ; for i in $(TARGETS); do \
	  cd $$i && \
	  $(MAKE) $(BUILDOPTS) BUILD=$(BUILD) VERBOSE=$(VERBOSE) && \
	  cd $$DIR; \
	done

install: dummy
	@$(MAKE) subs BUILDOPTS=install BUILD=$(BUILD) VERBOSE=1

depend: dummy
	@$(MAKE) subs BUILDOPTS=depend BUILD=$(BUILD)

test: dummy
	@$(MAKE) subs BUILDOPTS=test BUILD=$(BUILD) VERBOSE=1

testbuild: dummy
	@$(MAKE) subs BUILDOPTS=testbuild BUILD=$(BUILD) VERBOSE=1

relink: dummy
	@$(MAKE) subs BUILDOPTS=relink BUILD=$(BUILD)

clean: cleanlibs clean_targets

clean_targets:
	@-$(MAKE) subs BUILDOPTS=clean TARGETS="$(TARGETS)"

cleanlibs: dummy
	@-rm -f lib/lib*.a
	@-rm -f tmp/*.pch

cleanall: cleanlibs
	@-$(MAKE) subs BUILDOPTS=cleanall TARGETS="$(TARGETS)"

clean_source: cleanlibs
	@-$(MAKE) subs BUILDOPTS=clean_source TARGETS="$(TARGETS)"


dummy:

help:
	@cat $(ONETOOLROOT)/config/Makefile.help
