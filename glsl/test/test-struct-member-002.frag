#ifdef GL_ES
precision mediump float;
#endif


struct Ray {
	vec3 origin;
    vec3 direction;
};

struct Muda {
	vec3 morigin;
    vec3 mdirection;
};

struct Sphere {
	vec3 center;
	float radius;
};
	
void main( void ) {

    //Ray ray;
    Sphere s;

	gl_FragColor = vec4( s.center, s.radius );

}
