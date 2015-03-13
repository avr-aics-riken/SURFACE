#ifdef GL_ES
precision highp float;
#endif

uniform mat4 tmat;

void main(void)
{
    gl_FragColor = tmat[0];
}
