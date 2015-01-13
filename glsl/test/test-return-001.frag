#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;
varying vec3 mnormal;

void main(void) {
	
	int dep;
	raydepth(dep);
	if (dep > 8) {
		gl_FragColor = vec4(0,0,0,0);
		return; // Fail!!
	} else {
		vec3 p;
		vec3 n;
		vec3 dir;
		isectinfo(p, n, dir);
		
		vec3 pp = p + 0.01 * dir;
		
		vec4 col = vec4(0.0,0.0,0.0,0.0);
		float dist = trace(pp, dir, col);// NG
		//float dist = trace(pp, dir);// OK
		
		float y = 0.3;
		if (dist < -1000.0)
			y = 0.3 + col.r;
	
		gl_FragColor = vec4(y, y, y, 1);
	}
}
