#extension GL_OES_texture_3D : enable

#ifdef GL_ES
precision mediump float;
#endif

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main(void) {
    vec2 coord;
    coord.x = 0.0;
    coord.y = 0.0;
    gl_FragColor = vec4(resolution.x,0.0,0.0,1.0);
}
