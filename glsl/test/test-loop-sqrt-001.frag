#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main( void ) {
     const int MAXITER = 90;
     const int MAXITER_SQR= MAXITER*MAXITER;

     vec3 color = vec3(0);
     for (int i = 0; i < MAXITER; i++) {
	 color += float(MAXITER-i);
     }
     vec3 color3 = vec3(1.-1./sqrt(1.+color*(.09/float(MAXITER_SQR))));
     color3 *= color3;
     gl_FragColor = vec4(vec3(color3.r+color3.g+color3.b),1.);
}


