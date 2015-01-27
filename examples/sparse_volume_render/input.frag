#extension GL_OES_texture_3D : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler3D tex0;
uniform vec2      resolution;
uniform vec3      eye;
uniform vec3      lookat;
uniform vec3      up;

const float fov = 45.0;
const float num_steps = 64.0;

void
build_cameraframe(out vec3 corner, out vec3 u, out vec3 v)
{
    vec3 look = normalize(lookat - eye);
    vec3 uu = normalize(cross(look, up));
    vec3 vv = normalize(cross(look, uu));

    float flen = (0.5f * resolution.y / tan(0.5 * (fov * 3.14 / 180.0)));
    look = flen * look + eye;

    u = uu;
    v = vv;
    corner = look - 0.5 * (resolution.x * uu + resolution.y * vv);
}

// Simple transfer function
vec3 density_to_color(float dens)
{
    if (dens < 0.0) {
        return vec3(0.0, 0.0, 0.0);
    } else if (dens < 0.25) {
        float t = 4.0 * dens; // [0, 1]
        return vec3(t, 0.2*(1.0 - t), 0.0);
    } else if (dens < 0.5) {
        float t = 4.0 * (dens - 0.25); // [0, 1]
        return vec3(0.0, t, 1.0-t);
    } else if (dens < 0.75) {
        float t = 4.0 * (dens - 0.5); // [0, 1]
        return vec3(1.0-t, 0.0, t);
    } else {
        return vec3(0.0, 0.0, 0.0);
    }
}

int inside(vec3 p) {
    if ((p.x >= -1.0) && (p.x < 1.0) &&
        (p.y >= -1.0) && (p.y < 1.0) &&
        (p.z >= -1.0) && (p.z < 1.0)) {
        return 1;
    }
    return 0;
}

int
IntersectP(vec3 rayorg, vec3 raydir, vec3 pMin, vec3 pMax, out float hit0, out float hit1) {
    float t0 = -10000.0, t1 = 10000.0;
    hit0 = t0;
    hit1 = t1;

    vec3 tNear = (pMin - rayorg) / raydir;
    vec3 tFar  = (pMax - rayorg) / raydir;
    if (tNear.x > tFar.x) {
        float tmp = tNear.x;
        tNear.x = tFar.x;
        tFar.x = tmp;
    }
    t0 = max(tNear.x, t0);
    t1 = min(tFar.x, t1);

    if (tNear.y > tFar.y) {
        float tmp = tNear.y;
        tNear.y = tFar.y;
        tFar.y = tmp;
    }
    t0 = max(tNear.y, t0);
    t1 = min(tFar.y, t1);

    if (tNear.z > tFar.z) {
        float tmp = tNear.z;
        tNear.z = tFar.z;
        tFar.z = tmp;
    }
    t0 = max(tNear.z, t0);
    t1 = min(tFar.z, t1);

    if (t0 <= t1) {
        hit0 = t0;
        hit1 = t1;
        return 1;
    }
    return 0;
}

void  main(void) {

    vec2 coord = gl_FragCoord.xy / resolution;

    vec3 corner, u, v;
    build_cameraframe(corner, u, v);

    vec3 rayorg = eye;
    vec3 raydir = (corner + gl_FragCoord.x * u + gl_FragCoord.y * v) - eye;
    raydir = normalize(raydir);

    vec4 sum = vec4(0.0, 0.0, 0.0, 0.0);

    vec3 pmin = vec3(-1.0, -1.0, -1.0);
    vec3 pmax = vec3(1.0, 1.0, 1.0);
    float tmin, tmax;

    if (IntersectP(rayorg, raydir, pmin, pmax, tmin, tmax) < 1) { // no hit
        gl_FragColor = vec4(0.0, 0.0, 0.25, 1.0);
        return;
    }
 
    // raymarch.
    float t = tmin;
    float tstep = (tmax - tmin) / num_steps;
    float cnt = 0.0;
    while (cnt < float(num_steps)) {
        vec3 p = rayorg + t * raydir;
        // [-1, 1]^3
        if (inside(p) > 0) {
            // [-1, 1] - > [0, 1]^3
            vec4 dens = texture3D(tex0, 0.5 * p + 0.5);
            vec3 col = density_to_color(dens.x);
            sum += vec4(col, 1.0);
        }
        t += tstep;
        cnt += 1.0;
    }

    sum /= num_steps;

    gl_FragColor = vec4(sum.xyz, 1.0);
}
