#!/bin/sh

flags="--with-openmp --with-k-profile"
../tools/linux/premake4 --platform="k-cross" ${flags} gmake
