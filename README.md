# SURFACE

Scalable and Ubiquitous Rendering Framework for Advanced Computing Environments.

See `doc/INSTALL.md` for building/installing SURFACE.


## Build with CMake

### Supported platforms

* K(SPARC HPC-ACE)
* Linux x64
* MacOSX 10.8 or later
  * gcc 4.8 or later required for OpenMP build

### Requirements

* CMake 2.8 or later
* OpenMP-supported compiler(optional)
* MPI(optional)

### K cross compiling 

In the K frontend machine,

    $ /opt/local/bin/cmake -DCMAKE_INSTALL_PREFIX=<</path/to/install_dir>> -H. -DCMAKE_BUILD_TYPE=Release -Bbuild -DSURFACE_BUILD_K_CROSS_COMPILE=On
    $ make -C build
    $ make install

### MacOSX/Linux 

    $ cmake -DCMAKE_INSTALL_PREFIX=<</path/to/install_dir>> -H. -DCMAKE_BUILD_TYPE=Release -Bbuild
    $ make -C build
    $ make install


### Compile options

* SURFACE_BUILD_K_CROSS_COMPILE On/Off Set cross compile environment using Fujitsu cross compiler.
* SURFACE_BUILD_WITH_MPI On/Off Enable MPI.
* SURFACE_BUILD_WITH_OPENMP On/Off Enable OpenMP.

## Featrures

### Supported geometric primitives

* [x] Polygon(triangle)
* [x] Particle(rendered as sphere)
* [x] Tetrahedron
* [x] Line(rendered as tube)
* [ ] Curve primitives(ribbon, subdivision surface, etc)

