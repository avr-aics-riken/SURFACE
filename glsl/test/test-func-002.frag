#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

vec4 bora(vec2 a)
//vec4 bora()
{
    return vec4(a.x, 2, a.y, 4);
}

void main( void )
{
    gl_FragColor = bora(resolution);
}
