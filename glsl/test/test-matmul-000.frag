#ifdef LSGL_ES
#endif

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;
uniform mat4      mm;

void main( void )
{

    //mat2 ma = mat2(resolution.x, resolution.y, resolution.x * resolution.y, resolution.x + resolution.y);
    float a;
    vec4 b;
    vec4 c = mm * b;
    gl_FragColor = c;
}
