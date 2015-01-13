#ifdef GL_ES
precision mediump float;
#endif

varying vec4  myvar;

void main()
{
    gl_FragColor = myvar;
}
