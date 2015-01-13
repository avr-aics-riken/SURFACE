#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 resolution;

void main( void ) {
    
     const int MAXITER = 30;
     vec3 color = vec3(0);
     for (int i = 0; i < MAXITER; i++) {
	 float q = float(MAXITER-i);
   	 color += q /(color.x+0.1);
     }
     gl_FragColor = vec4(color,1.0);
}

