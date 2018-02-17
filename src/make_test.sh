#!/bin/bash
echo "  Starting make_test.sh"
echo `pwd`
echo

ARCH="x86_64-linux-gcc"
echo "  Configuring for $ARCH."
echo "ARCH=$ARCH" > config/arch
echo "COMP=gcc"   >> config/arch

mkdir targets
mkdir targets/$ARCH

echo
echo "  Building directory tree."
echo

for DIR in `find $SRCDIRS -type d -print`; do
  echo $DIR
  mkdir targets/$ARCH/$DIR
done

