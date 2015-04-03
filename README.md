# SURFACE

Scalable and Ubiquitous Rendering Framework for Advanced Computing Environments.

An software implementation of high performance ray tracer with GLES2.0 compatible API and GLSL shading language support.

## Requirements

* CMake 2.8 or later
* OpenMP-supported compiler(optional)
* MPI compiler(optional)

## Building with CMake

### Supported platforms

* [x] K/FX10(SPARC HPC-ACE)
  * Fujitsu C/C++ compiler
* [x] Linux x64
  * gcc 4.8 or later required for OpenMP build
  * Intel compiler
* [x] MacOSX 10.8 or later
  * gcc 4.8 or later required for OpenMP build
* [ ] Windows 7 x64
  * [ ] Visual Studio 2012
  * [x] Visual Studio 2013 Win64 (Library build only)

### K cross compiling 

In the K frontend node:

    $ /opt/local/bin/cmake -DCMAKE_INSTALL_PREFIX=dist -H. -DCMAKE_BUILD_TYPE=Release -Bbuild -DSURFACE_BUILD_K_CROSS_COMPILE=On
    $ make -C build
    $ make -C build install

### MacOSX/Linux 

    $ cmake -DCMAKE_INSTALL_PREFIX=dist -DBUILD_SHARED_LIBS=Off -H. -DCMAKE_BUILD_TYPE=Release -Bbuild
    $ make -C build
    $ make -C build install

### Windows 

    > vcbuild.bat

Solution file will be geneted in `build` directory.

### Compile options

* BUILD_SHARED_LIBS On/Off Build shared libs or static libs(default: shared)
* SURFACE_BUILD_K_CROSS_COMPILE On/Off Set cross compile environment using Fujitsu cross compiler.
* SURFACE_BUILD_WITH_MPI On/Off Enable MPI.
* SURFACE_BUILD_SCREEN_PARALLEL On/Off Build with screen parallel support?(Also need to enable SURFACE_BUILD_WITH_MPI option)
* SURFACE_BUILD_WITH_OPENMP On/Off Enable OpenMP.

## Setup GLSL compiler

SURFACE library calls GLSL compiler as an external program when `glShaderSource` and `glCompileShader` has been called. To make this work, Seting `GLSL_COMPILER` environment to point `glslc` binary is mandatory.

    $ export GLSL_COMPILER=/path/to/SURFACE/dist/glsl/glslc

(On Windows, the path should point to `glslc.bat`, not `glslc`)

Edit `glsl_config.py` to fit your environment(Installed at `/path/to/SURFACE/dist/glsl/glsl_config.py`)

Then run GLSL compiler to check it can produce shader module(`shader.so`)

    $ $GLSL_COMPILER /path/to/SURFACE/examples/tetra_render/input.frag 

## Featrures

* [x] Efficient ray tracing with GLES2.0 compatible API 
* [x] Support of GLES2.0 shading language(GLSL)
  * [ ] Vertex shader
  * [x] Fragment shader
* [x] Large scale data-parallel rendering using MPI
  * Up to 82,944 nodes on K computer confirmed.
* [ ] Sort-last image compositing using 234Compositor https://github.com/avr-aics-riken/234Compositor

### Supported geometric primitives

* [x] Polygon(triangle)
* [x] Particle(rendered as sphere)
* [x] Tetrahedron
* [x] Line(rendered as cylinder)
* [x] Volume
  * [x] Uniform volume
  * [x] Non-uniform volume
  * [x] Hierarchical volume(Sparse volume)

## Publications

* LSGL: Large Scale Graphics Library for Peta-Scale Computing Environments. HPG 2014 poster.
