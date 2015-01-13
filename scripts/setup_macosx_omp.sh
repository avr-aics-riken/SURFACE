#!/bin/sh

flags="--with-debugtrace --with-openmp"

../tools/macosx/premake4 ${flags} gmake
