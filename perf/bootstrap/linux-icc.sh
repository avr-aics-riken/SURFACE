#!/bin/sh

flags="--with-openmp"

../tools/linux64/premake4 ${flags} --platform=icc gmake
