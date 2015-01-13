#!/bin/sh

olddir=`pwd`

cd ~/src
diff -crN Mesa-9.0.1.org/src/glsl Mesa-9.0.1.lsgl/src/glsl > "${olddir}"/mesa_glsl_fix.patch

cd $olddir
