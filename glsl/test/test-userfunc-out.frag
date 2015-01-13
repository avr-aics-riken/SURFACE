#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec3      eyepos;
uniform vec3      eyedir;
uniform vec2      resolution;
varying vec3 mnormal;

void getV(out vec4 v)
{
	v = vec4(1.0, 1.0, 1.0, 1.0);
}


void main(void) {
	vec4 t;
	getV(t);
	gl_FragColor = t;
}
