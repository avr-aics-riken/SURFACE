#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

void main( void ) {

    vec3 eye;
    vec3 target;
    vec3 up;

    float fov = camerainfo(eye, target, up);

	gl_FragColor = vec4(eye, target.x);

}
