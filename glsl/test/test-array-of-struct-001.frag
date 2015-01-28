#ifdef GL_ES
precision mediump float;
#endif

uniform float t;

struct Sphere{
    float radius;
    vec3  position;
    vec3  color;
};

Sphere sphere[3];

void func(Sphere S, out float r){
    r = S.radius;
}

void main(void){

    // sphere init
    sphere[0].radius = 0.5;
    sphere[0].position = vec3(0.0, -0.5, sin(t));
    sphere[0].color = vec3(0.2, 0.9, 0.7);
    sphere[1].radius = 1.0;
    sphere[1].position = vec3(2.0, 0.0, cos(t * 0.666));
    sphere[1].color = vec3(0.8, 0.9, 0.9);
    sphere[2].radius = 1.5;
    sphere[2].position = vec3(-2.0, 0.5, cos(t * 0.333));
    sphere[2].color = vec3(0.0, 0.0, 1.0);

    float r0;
    func(sphere[0], r0);
    float r1;
    func(sphere[1], r1);
    float r2;
    func(sphere[2], r2);

    gl_FragColor = vec4(r0, r1, r2, 1.0);
}
