#extension GL_LSGL_random : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 resolution;

uniform float aa;

void main()
{
    float b;
    random(b);

    gl_FragColor = vec4(b); 
}
