# SURFACE

Scalable and Ubiquitous Rendering Framework for Advanced Computing Environments.

See `doc/INSTALL.md` for building/installing SURFACE.


## Build with CMake

### Requirements

* CMake 2.8 or later

### Build

    $ cmake -H. -DCMAKE_BUILD_TYPE=Release -Bbuild
    $ make -C build

## Featrures

### Supported geometric primitives

* [x] Polygon(triangle)
* [x] Particle(rendered as sphere)
* [x] Tetrahedron
* [x] Line(rendered as tube)
* [ ] Curve primitives(ribbon, subdivision surface, etc)

