#!/bin/sh

flags="--with-debugtrace --with-mpi --with-screen-parallel"

../tools/macosx/premake4 ${flags} gmake

# change CC, CXX compiler to mpicc/mpicxx
# MacOS's sed doesn't handle newline well, so use GNU sed.
gsed -i -e "s/ifndef CC/CC=mpicc\nifndef CC/g" LSGLES.make
gsed -i -e 's/ifndef CXX/CXX=mpicxx\nifndef CXX/g' LSGLES.make
