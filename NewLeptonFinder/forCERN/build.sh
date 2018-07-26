#!/bin/sh

if [ -d "$HOME/NewLeptonFinder/forCERN/build" ]; then
   rm -rf "$HOME/NewLeptonFinder/forCERN/build"
fi

mkdir "$HOME/NewLeptonFinder/forCERN/build"
cd "$HOME/NewLeptonFinder/forCERN/build"
cmake -C $ILCSOFT/ILCSoft.cmake ..
make
make install
