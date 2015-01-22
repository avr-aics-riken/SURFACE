#!/bin/sh

olddir=`pwd`
rm -rf build
mkdir build
cd build
tar -jxvf ../../deps/MesaLib-9.0.1.tar.bz2
cd Mesa-9.0.1
patch -p1 < ../../mesa_glsl_fix.patch
cd $olddir
