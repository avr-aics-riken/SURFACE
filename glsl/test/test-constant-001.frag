#extension GL_OES_texture_3D : enable

#ifdef GL_ES
precision mediump float;
#endif


//varying vec3  Normal;
//varying vec3  LightDir;
//varying vec3  Eye;

//varying vec3  Reflect;
//varying float Ratio;

uniform sampler3D tex0;
//uniform samplerCube Cubemap;
//uniform float myval[4];
//uniform vec3 myval02[4];
//uniform bvec4 myvec[4];
uniform vec2 resolution;
//uniform vec2 coord;

uniform float aa;

void main()
{
    //gl_FragColor.xyz = Normal;
    //float t = 0.0;

    vec4 uPos;
    uPos.x = 1.0;
    uPos.y = 0.0;
    uPos.z = 0.0;
    uPos.w = 1.0;
    gl_FragColor = uPos;

    //float t = gl_FragCoord.y;
    //gl_FragColor.x = 0.1f * t;
    //if (gl_FragCoord.y < 0.5) {
    //    gl_FragColor.x = gl_FragCoord.z;
    //}

    //vec3 coord = gl_FragCoord.xyz;
    //gl_FragColor = texture3D(tex0, coord);

    //float yy = gl_FragCoord.y / resolution.y;

    //while (t < 1.0) {
    //    gl_FragColor.xyz = 0.1f * gl_FragCoord.xyz;
        //t += cos(gl_FragCoord.x);
    //}
    //vec3 In = normalize(Eye);
    //gl_FragColor.xyz = In;
    //gl_FragColor.x = aa;
    //gl_FragColor.x = yy;
    //gl_FragColor.x = gl_FragCoord.x;
    //gl_FragColor.xyz = cos(gl_FragCoord.xyz);
    //gl_FragColor.xyz = m;
}
