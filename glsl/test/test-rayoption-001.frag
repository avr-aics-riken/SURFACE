#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 resolution;

uniform float aa;

void main()
{
    int prev_double_sided;
    int double_sided = 0;

    rayoption(prev_double_sided, double_sided);

    gl_FragColor = vec4(float(prev_double_sided)); 
}
