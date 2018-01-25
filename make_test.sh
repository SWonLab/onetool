#!/bin/bash
set -ev
pwd=`${pwd}`
echo pwd
if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then
  cd src
  echo "in src directory"
  pwd=`${pwd}`
  echo pwd
  #./configure
  #cd c++
  #make
fi
