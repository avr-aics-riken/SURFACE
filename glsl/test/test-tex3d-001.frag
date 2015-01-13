#extension GL_OES_texture_3D : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler3D tex0;

void main()
{
    vec3 coord = gl_FragCoord.xyz;
    gl_FragColor = texture3D(tex0, coord);
}
