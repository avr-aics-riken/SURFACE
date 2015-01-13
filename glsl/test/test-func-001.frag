#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

float testfunc(vec2 p)
{
    return 1.0;
}

void main(void) {
    vec2 coord;
    coord = gl_FragCoord.xy / resolution;

    vec2 p;
    p.x = testfunc(coord);
    p.y = testfunc(coord);
    gl_FragColor = vec4(p.x,p.y,1.0,1.0);
}
