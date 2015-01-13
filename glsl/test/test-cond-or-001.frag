#ifdef GL_ES
precision mediump float;
#endif

void main()
{
    if ((gl_FragCoord.x > 0.5) || (gl_FragCoord.y > 0.5)) {
        gl_FragColor.x = 1.0;
    } else {
        gl_FragColor.x = 0.1;
    }
}
