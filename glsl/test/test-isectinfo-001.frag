#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

void main( void ) {

    vec3 org;
    vec3 norm;
    vec3 dir;

    isectinfo(org, norm, dir);

	gl_FragColor = vec4(org, norm.x);

}
