#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec3      eyepos;
uniform vec3      eyedir;
uniform vec2      resolution;
varying vec3 mnormal;


void main(void) {
	float e = dot(normalize(mnormal),vec3(0,0,1)); // correct 0.0-1.0
#if 1
	e = max(0.0,e); // error! All to 0.0
#else
	// to expect same result with e = max(0.0,e)
	if (e < 0.0)
		e = 0.0;
#endif
	gl_FragColor = vec4(e,e,e,1.0);

}
