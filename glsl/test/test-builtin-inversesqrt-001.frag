#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main(void) {
    vec2 coord;
    coord = gl_FragCoord.xy / resolution;

    vec2 p;
    p.x = inversesqrt(coord.x);
    gl_FragColor = vec4(p.x,p.y,1.0,1.0);
}
