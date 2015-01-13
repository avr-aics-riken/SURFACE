#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;
//uniform vec4      vval;
//varying vec2 texcoord;
varying vec3 mnormal;

vec3 degamma(vec3 col)
{
	return pow(col,2.2); //NG
	//return vec3(pow(col.r,2.2),pow(col.g,2.2),pow(col.b,2.2));//OK
}

void main(void) {
	vec3 N = normalize(mnormal);
	float d = (N.y*.5+.5)*.8+.2;
    gl_FragColor = vec4(degamma(vec3(d,d,d)), 1.0);
}
