# LSGL extensions

* LSGL_trace
* LSGL_random

## LSGL_trace

```
    int numIntersects(out int num);
```

Returns the number of intersection points. For polygonal object, this is usually 1, and 1 or 2 for tetrahedron primitive.

```
    float queryIntersect();
```

## LSGL_random

```
    float random(out float random);
```
