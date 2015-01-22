# non-uniform volume rendering.

We represent non-uniform volume rendering using `lsglTexCoordRemap` custom function to add texture coordinate remapping feature when looking up textel in GLSL.
`lsglTexCoordRemap` is currently only valid for 3D textures(3D volumes).

## Limitation.

Sizeo of remap table should be equal or larger than volume resolution(e.g. if your volume has 64 voxels in Y direction, the size of remap table should be 64 or larger).

Coordinate remapping is specified with 1D remap table for each direction, but we believe this is enough for doing non-uniform volume rendering.
