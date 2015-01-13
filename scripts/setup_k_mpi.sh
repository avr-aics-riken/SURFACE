#!/bin/sh

flags="--with-mpi --with-openmp"
../tools/linux/premake4 --platform=k-cross-mpi ${flags} gmake

# change CC, CXX compiler to fccpx, FCCpx
#sed -i -e 's/ifndef CC/CC=mpifccpx\nifndef CC/g' LSGLES.make
#sed -i -e 's/ifndef CXX/CXX=mpiFCCpx\nifndef CXX/g' LSGLES.make
#sed -i -e 's/ifndef CC/CC=mpifccpx\nifndef CC/g' LSGLESTest.make
#sed -i -e 's/ifndef CXX/CXX=mpiFCCpx\nifndef CXX/g' LSGLESTest.make

## change include file
#sed -i -e 's/-f LSGLES.make/-f LSGLES.make.k/g' Makefile
#cp Makefile Makefile.k
#cp LSGLES.make LSGLES.make.k
