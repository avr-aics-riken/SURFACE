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

vec3 myfunc(vec3 org)
{
    return org + org;
}
	
void main( void ) {

    //Ray ray;
    Sphere s;
    vec3 aa = myfunc(s.center);

	gl_FragColor = vec4( aa, 1.0 );

}
