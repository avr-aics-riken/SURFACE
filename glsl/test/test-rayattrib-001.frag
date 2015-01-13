#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 resolution;

uniform float aa;

void main()
{
    float a;
    a = rayattrib(a);
    gl_FragColor = vec4(a); 
}
