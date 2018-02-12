#!/bin/bash
set -ev
pwdnow=`${pwd}`
echo pwdnow
if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then
  cd src
  echo "in src directory"
  pwdsrc=`${pwd}`
  echo pwdsrc
  #./configure
  #cd c++
  #make
fi
