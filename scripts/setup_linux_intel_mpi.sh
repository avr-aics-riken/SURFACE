#!/bin/sh

../tools/linux/premake4 --with-mpi --with-openmp gmake

# change CC, CXX compiler
sed -i -e 's/ifndef CC/CC=mpifcc\nifndef CC/g' LSGLES.make
sed -i -e 's/ifndef CXX/CXX=mpicxx\nifndef CXX/g' LSGLES.make
sed -i -e 's/ifndef CC/CC=mpifcc\nifndef CC/g' LSGLESTest.make
sed -i -e 's/ifndef CXX/CXX=mpicxx\nifndef CXX/g' LSGLESTest.make
