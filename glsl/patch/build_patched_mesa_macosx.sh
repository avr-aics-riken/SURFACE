#!/bin/sh

olddir=`pwd`
cd build/Mesa-9.0.1/src/glsl/
scons -u -c
scons -u debug=0 build=release
if [ -f ${olddir}/build/Mesa-9.0.1/build/darwin-x86_64/glsl/glsl2 ]; then
	cp ${olddir}/build/Mesa-9.0.1/build/darwin-x86_64/glsl/glsl2 ${olddir}/../bin/macosx64/glsl_compiler
fi
cd $olddir
