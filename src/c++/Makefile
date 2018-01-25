#==========================================================================
# Makefile for oneTool
#
#  Author:  Yeunjoo E. Song (yeunjoo.song@gmail.com)
#
#  History: Initial implementation                        - yes Nov 24 2014
#
#  Notes:   oneTool suite-wide make
#
#  Copyright (c) 2014  Sungho Won
#    All Rights Reserved
#==========================================================================

VERSION    = 2
PATCHLEVEL = 0
SUBLEVEL   = 0
VERBOSE    = 3
BUILDOPTS  =
CXXFLAGS   = -g

ONETOOLROOT=..
include $(ONETOOLROOT)/config/Global.make

CORE        = global utils
IO          = input output
ANALYSES    = analyses sage merlin
APPLICATION = app

TARGETS = $(CORE) $(IO) $(ANALYSES) $(APPLICATION)

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

docs: dummy
	@echo "[lib categories]"                 >  docs/~~docdict
	@echo "CORE:          $(CORE)"           >> docs/~~docdict
	@echo "DATA:          $(DATA)"           >> docs/~~docdict
	@echo "IO:            $(IO)"             >> docs/~~docdict
	@echo "INFO:          $(INFO)"           >> docs/~~docdict
	@echo "ANALYSIS:      $(ANALYSIS)"       >> docs/~~docdict
	@echo "APPLICATION:   $(APPLICATION)"    >> docs/~~docdict

	@$(MAKE) subs BUILDOPTS=docs BUILD=$(BUILD) VERBOSE=1

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
