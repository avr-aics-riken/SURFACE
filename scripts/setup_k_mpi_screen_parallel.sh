#!/bin/sh

flags="--with-debugtrace --with-screen-parallel --with-mpi"
../tools/linux/premake4 --platform="k-cross-mpi" ${flags} gmake
