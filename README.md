# SURFACE

Scalable and Ubiquitous Rendering Framework for Advanced Computing Environments.

## Build with CMake

### Supported platforms

* [x] K(SPARC HPC-ACE)
* [x] Linux x64
* [x] MacOSX 10.8 or later
  * gcc 4.8 or later required for OpenMP build
* [ ] Windows 7 x64

### Requirements

* CMake 2.8 or later
* OpenMP-supported compiler(optional)
* MPI(optional)

### K cross compiling 

In the K frontend machine,

    $ /opt/local/bin/cmake -DCMAKE_INSTALL_PREFIX=dist -H. -DCMAKE_BUILD_TYPE=Release -Bbuild -DSURFACE_BUILD_K_CROSS_COMPILE=On
    $ make -C build
    $ make -C build install

### MacOSX/Linux 

    $ cmake -DCMAKE_INSTALL_PREFIX=dist -DBUILD_SHARED_LIBS=Off -H. -DCMAKE_BUILD_TYPE=Release -Bbuild
    $ make -C build
    $ make -C build install

### Compile options

* BUILD_SHARED_LIBS On/Off Build shared libs or static libs(default: shared)
* SURFACE_BUILD_K_CROSS_COMPILE On/Off Set cross compile environment using Fujitsu cross compiler.
* SURFACE_BUILD_WITH_MPI On/Off Enable MPI.
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
* [x] Large scale data-parallel rendering using MPI
  * Up to 82,944 nodes on K confirmed.

### Supported geometric primitives

* [x] Polygon(triangle)
* [x] Particle(rendered as sphere)
* [x] Tetrahedron
  * [ ] Remove rendering artifact
* [x] Line(rendered as tube)
* [x] Volume
  * [ ] Uniform volume
  * [ ] Non-uniform volume
  * [ ] Hierarchical volume
