#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
varying vec2      texcoord;

void main()
{
    gl_FragColor = texture2D(tex0, texcoord);
}
