#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main( void )
{

    float time = 0.5;
    vec2 pos = (gl_FragCoord.xy*2.0 - resolution) / resolution.y;
    vec3 p = vec3(pos,1.0);
    vec3 n;
    n = floor(p - vec3(0.0,0.0,0.5));
    gl_FragColor = vec4(n,1.0);
}