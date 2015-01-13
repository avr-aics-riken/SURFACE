#!/bin/sh

../tools/linux/premake4 --with-openmp gmake

# change CC, CXX compiler to icc, icpc
sed -i -e 's/ifndef CC/CC=icc\nifndef CC/g' LSGLES.make
sed -i -e 's/ifndef CXX/CXX=icpc\nifndef CXX/g' LSGLES.make
sed -i -e 's/ifndef CC/CC=icc\nifndef CC/g' LSGLESTest.make
sed -i -e 's/ifndef CXX/CXX=icpc\nifndef CXX/g' LSGLESTest.make
