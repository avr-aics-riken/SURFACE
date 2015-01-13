#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main( void )
{

    float time = 0.5;
    vec2 pos = (gl_FragCoord.xy*2.0 - resolution) / resolution.y;
    vec3 n;
    n.x = floor(pos.x);
    n.y = floor(pos.y);
    n.z = floor(-0.5 + pos.x);
    gl_FragColor = vec4(n,1.0);
}