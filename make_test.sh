#!/bin/bash
set -ev
if [ "${TRAVIS_PULL_REQUEST}" = "false" ]; then
  cd src
  ./configure
  cd c++
  make
fi
