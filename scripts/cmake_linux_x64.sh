#!/bin/sh

if [ -z "${CMAKE_BIN+x}" ]; then
  CMAKE_BIN=cmake28
fi

echo ${CMAKE_BIN}
${CMAKE_BIN} -H. -DCMAKE_INSTALL_PREFIX=dist -DBUILD_SHARED_LIBS=On -DCMAKE_BUILD_TYPE=Release -Bbuild
