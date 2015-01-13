#ifdef LSGL_ES
float MAX_ITER = 100.0;
#endif

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main( void )
{

    float time = 0.5;
    float i;
    discard;
    for (i = 0.0; i < MAX_ITER; i++) {
        gl_FragColor.x += 0.1;
    }
}
