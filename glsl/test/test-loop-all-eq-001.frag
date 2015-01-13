#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main( void ) {
     vec3 color = vec3(0);
     const int MAXITER = 30;
     for (int i = 0; i != MAXITER; i++) {
         color += 0.1*float(i/MAXITER);
     }
     gl_FragColor = vec4(color,1.0);
}
										
										


