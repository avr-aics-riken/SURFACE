#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 resolution;

uniform float aa;

void main()
{
    vec3 p, n;
    vec4 col;
    float myattrib;
    float a;
    a = trace(p, n, col, myattrib);
    gl_FragColor = vec4(a); 
}
