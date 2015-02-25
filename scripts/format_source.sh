#!/bin/sh

curdir=`pwd`

cd gles;
clang-format -i gles_*.cc
clang-format -i gles_*.h
clang-format -i ../render/*.cc
clang-format -i ../render/*.h

cd ${curdir}
