#!/bin/bash
echo
echo "  Starting configure & make process..."
echo `pwd`
echo

ONETOOLROOT=.
echo $ONETOOLROOT

SRCDIRS="c++"
echo $SRCDIRS

hardware=`uname -m`
echo $hardware

osname=`uname -s`
echo $osname

release=`uname -r`
echo $release

ARCH="x86_64-linux-gcc"
echo
echo "  Configuring for $ARCH..."
echo "ARCH=$ARCH" > config/arch
echo "COMP=gcc"   >> config/arch

mkdir targets
mkdir targets/$ARCH

echo
echo "  Building directory tree..."
echo

for DIR in `find $SRCDIRS -type d -print`; do
  #echo $DIR
  mkdir targets/$ARCH/$DIR
done

echo
echo "  Making ONETOOL..."
echo

cd $SRCDIRS
echo `pwd`
make
cd ..

echo
echo "  Done."
echo
