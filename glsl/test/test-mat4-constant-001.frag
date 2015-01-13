#ifdef LSGL_ES
#endif

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

mat4 transmat(mat4 m)
{
        return mat4(1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);
}

void main( void )
{
    mat4 m0, m1;
    m0 = transmat(m1);

    float time = 0.5;
    float i;
    float MAX_ITER = 100.0;
    discard;
    for (i = 0.0; i < MAX_ITER; i++) {
        gl_FragColor.x += 0.1;
    }
    gl_FragColor.y = m0[0][1];
}
