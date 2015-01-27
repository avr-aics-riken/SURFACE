#ifdef LSGL_ES
#endif

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main( void )
{

    mat2 ma = mat2(resolution.x, resolution.y, resolution.x * resolution.y, resolution.x + resolution.y);

    float time = 0.5;
    float i;
    float MAX_ITER = 100.0;
    discard;
    for (i = 0.0; i < MAX_ITER; i++) {
        gl_FragColor.x += 0.1 + ma[0][0];
    }
}
