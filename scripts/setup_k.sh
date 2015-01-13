#!/bin/sh

flags="--with-openmp"
../tools/linux/premake4 --platform="k-cross" ${flags} gmake
