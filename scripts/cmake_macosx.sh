#!/bin/sh

if [ -z "${CMAKE_BIN+x}" ]; then
  CMAKE_BIN=cmake
fi

if [ -z "${CXX+x}" ]; then
  CXX=g++
fi

if [ -z "${CC+x}" ]; then
  CC=gcc
fi

${CMAKE_BIN} -DCMAKE_CXX_COMPILER=${CXX} -DCMAKE_C_COMPILER=${CC} -DCMAKE_INSTALL_PREFIX=dist -H. -DBUILD_SHARED_LIBS=Off -DCMAKE_BUILD_TYPE=Release -Bbuild
