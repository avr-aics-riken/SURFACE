#ifdef GL_ES
precision mediump float;
#endif

varying vec3  Normal;
varying vec3  LightDir;
varying vec3  Eye;

//varying vec3  Reflect;
//varying float Ratio;

//uniform samplerCube Cubemap;
//uniform float myval[4];
//uniform vec3 myval02[4];
//uniform bvec4 myvec[4];

void main()
{
    if (gl_FragCoord.x > 0.5) {
        gl_FragColor.x = 1.0;
    } else {
        gl_FragColor.x = 0.1;
    }
}
