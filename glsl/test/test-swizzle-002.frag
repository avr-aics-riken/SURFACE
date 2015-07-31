#ifdef GL_ES
precision highp float;
#endif

uniform vec2 iResolution;
uniform float iGlobalTime;

vec2 rotate(vec2 p, float a)
{
    return vec2(p.x * cos(a) - p.y * sin(a), p.x * sin(a) + p.y * cos(a));
}

void main(void)
{
    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    
    uv = 2.0 * uv - 1.0;
    uv.x *= iResolution.x / iResolution.y;
    
    float camtm = iGlobalTime * 0.15;
    vec3 ro = vec3(cos(camtm), 0.0, camtm);
    vec3 rd = normalize(vec3(uv, 1.2));
    rd.xz = rotate(rd.xz, sin(camtm) * 0.4);
    rd.yz = rotate(rd.yz, sin(camtm * 1.3) * 0.4);
    
    gl_FragColor = vec4(rd.xyz, 1.0);
}
