#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

void main( void ) {

    vec3 org;
    vec3 dir;
    vec4 col;
	gl_FragColor = vec4(trace(org, dir, col), 1.0, 1.0, 1.0);

}
