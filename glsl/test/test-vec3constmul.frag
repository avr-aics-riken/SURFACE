#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;
/*
float rand(vec3 n, float res)
{
  n = floor(n*res);
  return fract(sin((n.x+n.y*1e2+n.z*1e4)*1e-4)*1e5);
}

float map( vec3 p )
{
    float d = rand(floor(p-vec3(0.,0.,0.5)), 1.);
    p = mod(p,vec3(1.0))-0.5;
    return length(p.xy)-(d*.2);
}
*/
void main( void )
{
    float time = 0.5;
    vec2 pos = (gl_FragCoord.xy*2.0 - resolution) / resolution.y;
    vec3 camTarget = vec3(sin(time), 0.0, .0);

    gl_FragColor = vec4(camTarget * pos.x,1.0);
}