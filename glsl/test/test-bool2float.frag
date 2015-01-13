#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main( void ) {
     vec2 p = -1.0 + 2.0 * gl_FragCoord.xy / resolution;
     gl_FragColor = vec4(p.x > 0.0,0,0,1.0);
}


