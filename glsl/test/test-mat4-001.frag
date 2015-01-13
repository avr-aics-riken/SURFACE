#ifdef LSGL_ES
#endif

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

mat4 transmat(mat4 m)
{
        return mat4(vec4(m[0][0],m[1][0],m[2][0],m[3][0]),
                                vec4(m[0][1],m[1][1],m[2][1],m[3][1]),
                                vec4(m[0][2],m[1][2],m[2][2],m[3][2]),
                                vec4(m[0][3],m[1][3],m[2][3],m[3][3]));
}

void main( void )
{

    float time = 0.5;
    float i;
    float MAX_ITER = 100.0;
    discard;
    for (i = 0.0; i < MAX_ITER; i++) {
        gl_FragColor.x += 0.1;
    }
}
