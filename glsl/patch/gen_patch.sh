#!/bin/sh

olddir=`pwd`
cd build
rm -rf Mesa-9.0.1.org
mkdir Mesa-9.0.1.org
tar -jxvf ../../deps/MesaLib-9.0.1.tar.bz2 -C Mesa-9.0.1.org --strip-components=1
diff -crN Mesa-9.0.1.org/src/glsl Mesa-9.0.1/src/glsl > "${olddir}"/mesa_glsl_fix.patch
cd $olddir
