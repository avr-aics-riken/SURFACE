#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 resolution;

vec3 field(vec3 p) {
     p = abs(fract(p)-.5);
     p *= p;
     return sqrt(p+p.yzx*p.zzy)-.015;
}

void main( void ) {
     vec3 color = vec3(0.0);
     gl_FragColor = vec4(color,1.0);
}