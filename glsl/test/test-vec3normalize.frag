#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;

void main( void )
{
    float time = 0.5;
    vec2 pos = (gl_FragCoord.xy*2.0 - resolution) / resolution.y;

    vec3 camPos = vec3(cos(time*0.3), sin(time*0.3), 3.5);
    vec3 camTarget = vec3(0.0, 0.0, 0.0);

    vec3 camDir = normalize(camTarget-camPos);
    vec3 camUp  = normalize(vec3(0.0, 1.0, 0.0));
    vec3 camSide = cross(camDir, camUp);
    float focus = 1.8;

   vec3 rayDir;
   rayDir.x = camSide.x * pos.x + camUp.x * pos.y + camDir.x * focus;
   rayDir.y = camSide.y * pos.x + camUp.y * pos.y + camDir.y * focus;
   rayDir.z = camSide.z * pos.x + camUp.z * pos.y + camDir.z * focus;
   rayDir = normalize(rayDir);

   gl_FragColor = vec4(rayDir,1.0);
}