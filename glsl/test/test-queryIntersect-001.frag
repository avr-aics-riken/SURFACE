#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

void main( void ) {

    vec3 org;
    vec3 Ng, Ns;
    vec3 Tn, Bn;
    vec3 dir;
    vec2 bary;
    float t = queryIntersect(0, org, Ng, Ns, Tn, Bn, dir, bary);
	gl_FragColor = vec4(float(t), 1.0, 1.0, 1.0);

}
