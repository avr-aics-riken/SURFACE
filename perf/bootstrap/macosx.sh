#!/bin/sh

flags="--with-openmp"

../tools/macosx/premake4 ${flags} gmake

# change CC, CXX compiler to clang, clang++
gsed -i -e 's/ifndef CC/CC=gcc-4.8\nifndef CC/g' RenderBenchmark.make
gsed -i -e 's/ifndef CXX/CXX=g++-4.8\nifndef CXX/g' RenderBenchmark.make

# change CC, CXX compiler to clang, clang++
gsed -i -e 's/ifndef CC/CC=gcc-4.8\nifndef CC/g' LibRender.make
gsed -i -e 's/ifndef CXX/CXX=g++-4.8\nifndef CXX/g' LibRender.make
