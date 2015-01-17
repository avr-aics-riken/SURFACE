PROJECT(SURFACE)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

SET(SURFACE_VERSION 0.9.0)
SET(SURFACE_SOVERSION 0)   # shared lib

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the build type: Debug Release RelWithDebInfo MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(SURFACE_BINARY_DIR ${PROJECT_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${SURFACE_BINARY_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${SURFACE_BINARY_DIR})
LINK_DIRECTORIES(${LIBRARY_OUTPUT_PATH})

MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

# Build options
SET(SURFACE_BUILD_UNIT_TEST OFF CACHE BOOL "Build unit test?")
SET(SURFACE_BUILD_K_CROSS_COMPILE OFF CACHE BOOL "Build for K/FX10 target with corss compiling?")
SET(SURFACE_BUILD_WITH_MPI OFF CACHE BOOL "Build with MPI support?")

# K/FX10 native environment
IF (CMAKE_SYSTEM_PROCESSOR MATCHES "s64fx")

  # Use native Fujitsu compiler
  SET(CMAKE_C_COMPILER fcc)
  SET(CMAKE_CXX_COMPILER FCC)

  # Prevent CMake from adding GNU-specific linker flags (-rdynamic)
  SET(CMAKE_C_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu C cross-compiler" FORCE)
  SET(CMAKE_CXX_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu C++ cross-compiler" FORCE)

  # FindOpenMP.cmake does not try -Kopenmp,but the package will try specific
  # flags based on the compier ID.
  SET(OMP_FLAG_Fujitsu "-Kopenmp")
ENDIF()

# K cross compiling
IF (SURFACE_BUILD_K_CROSS_COMPILE)

  # Use Fujitsu cross compiler
  SET(CMAKE_C_COMPILER fccpx)
  SET(CMAKE_CXX_COMPILER FCCpx)

  # Prevent CMake from adding GNU-specific linker flags (-rdynamic)
  SET(CMAKE_C_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu C cross-compiler" FORCE)
  SET(CMAKE_CXX_COMPILER_ID "Fujitsu" CACHE STRING "Fujitsu C++ cross-compiler" FORCE)

  # gcc compatible flag(-Xg)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xg")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xg")

  # OpenMP flag(-Kopenmp)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kopenmp")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kopenmp")

  IF (CMAKE_BUILD_TYPE MATCHES "Release")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast")
  ENDIF()

  # FindOpenMP.cmake does not try -Kopenmp,but the package will try specific
  # flags based on the compier ID.
  SET(OMP_FLAG_Fujitsu "-Kopenmp")
  
ENDIF()

# Modules
ADD_SUBDIRECTORY(gles gles)
#ADD_SUBDIRECTORY(examples examples)
