#!/bin/bash
echo
echo "  Starting configure & make process..."
echo `pwd`
echo

ONETOOLROOT="src"
cd $ONETOOLROOT
echo `pwd`
./configure

SRCDIRS="c++"
cd $SRCDIRS
echo `pwd`
make
cd ../..

echo
echo "  Done."
echo
