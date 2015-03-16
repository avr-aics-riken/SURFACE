#ifdef GL_ES
precision highp float;
#endif

void main() {
    vec3 I, N;
    vec3 R = refract(I, N, 1.22);
}

