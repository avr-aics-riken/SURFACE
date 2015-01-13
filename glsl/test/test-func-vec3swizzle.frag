#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;
varying vec3 mnormal;

vec3 map(vec3 dir)
{
	float r = dir.y;
	float l = max(0.0,r);
	return vec3(l, l, l);
}

void main(void) {
	//float r = mnormal.y;
	//float l = max(0.0,r);
	//vec3 l3= vec3(l, l, l); // <- no error
	vec3 l3 = map(mnormal);
    gl_FragColor = vec4(l3, 1.0);
}
