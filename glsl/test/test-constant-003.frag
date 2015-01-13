#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main(void) {
    vec2 coord;
    coord = gl_FragCoord.xy / resolution;

    //gl_FragColor = texture2D(tex0,vec2(coord.x + 0.1, coord.y + 0.1)); // Successed
    //gl_FragColor = texture2D(tex0,vec2(-0.1 + coord.x, -0.1 + coord.y)); // calc miss
    gl_FragColor = texture2D(tex0,vec2(coord.x - 0.1, coord.y - 0.1)); // Compiler Fail
}
