#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 resolution;

const float intensity = 0.5;

void main()
{
    vec2 c = vec2(gl_FragCoord.xy/resolution);
    c*=intensity;
    gl_FragColor = vec4(c, 0.0, 1.0);
}
