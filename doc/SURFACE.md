# SURFACE

T.B.W.

## SURFACE extensions

### GLES extensions

#### Point primitive

* `uniform float lsgl_PointSize` specify size(width) for all point primitives.
* `varying float lsgl_PointSize` specify size(width) for each point primitive.

#### Line primitive

* glLineWidth() to specify width for all line primitives.
* `uniform int lslg_LineCap` specify whether cap line(cylinder) primitive or not. Default value is 1(true).
* `varying float lsgl_LineWidth` specify size(width) for each line primitive.


### GLSL extensions

Since SURFACE doesn't support vertex shader at this time, predefined uniform variables are provided for vertex transform(like `GL_MODEVIEW` matrix in OpenGL).

### Uniforms

```
uniform mat4 lsgl_World; // object to world matrix
``` 

SURFACE use these matrix to transform object and set camera view


And there's predefined uniform variables which could be accessible from the fragment shader.

```
uniform mat4 lsgl_World; // object to world matrix
uniform mat4 lsgl_WorldInverse; // inverse of object to world matrix
uniform mat4 lsgl_WorldInverseTranspose; // inverse transpose of object to world matrix
``` 

#### Built-in functions

```
int numIntersections(out int n);
```

Retruns the number of intersections

T.B.W.
