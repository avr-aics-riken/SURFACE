#extension GL_LSGL_trace : enable

#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D tex0;
uniform vec2      resolution;
//uniform vec4      vval;
//varying vec2 texcoord;
//varying vec3  mnormal;
varying float matID;
varying vec3 colors;

float frac(float x) { return mod(x,1.0); }
float noise(vec2 x, float t) { return frac(sin(dot(x,vec2(1.38984*sin(t),1.13233*cos(t))))*653758.5453); }

#define PI 3.1415926

vec3 matToColor(float m)
{
    // NOTE: Replace following code with 1D texture fetch would be better idea!

    if (m < 0.0) { // Unsupported material color
        return vec3(0.75, 0.75, 1.0);
    }

    if (m < 1.1) { // Hydrogen
        return vec3(1.0, 1.0, 1.0);
    } else if (m < 2.1) { // Hellium
        return vec3(0.85, 1.0, 1.0);
    } else if (m < 3.1) { // Lithium
        return vec3(0.8, 0.5, 1.0);
    } else if (m < 4.1) { // Beryllium
        return vec3(0.76, 1.0, 0.0);
    } else if (m < 5.1) { // Boron
        return vec3(1.0, 0.7, 0.7);
    } else if (m < 6.1) { // Carbon
        return vec3(0.56, 0.56, 0.56);
    } else if (m < 7.1) { // Nitrogen
        return vec3(0.18, 0.31, 0.97);
    } else if (m < 8.1) { // Oxygen
        return vec3(1.0, 0.05, 0.05);
    } else if (m < 9.1) { // Fluorine
        return vec3(0.56, 0.87, 0.31);
    } else if (m > 9999.0) { // H2O
        return vec3(0.86, 0.2, 0.2);
    }

    // todo
    return vec3(0.5, 0.5, 0.5);
}


void orthoBasis(out vec3 basis0,out vec3 basis1,out vec3 basis2, vec3 n)
{
	basis2 = vec3(n.x, n.y, n.z);
	basis1 = vec3(0.0, 0.0, 0.0);

	if ((n.x < 0.6) && (n.x > -0.6))
	basis1.x = 1.0;
	else if ((n.y < 0.6) && (n.y > -0.6))
	basis1.y = 1.0;
	else if ((n.z < 0.6) && (n.z > -0.6))
	basis1.z = 1.0;
	else
	basis1.x = 1.0;

	basis0 = cross(basis1, basis2);
	basis0 = normalize(basis0);

	basis1 = cross(basis2, basis0);
	basis1 = normalize(basis1);
}

void main(void) {

	
    vec3 lightdir = normalize(vec3(1, 1, -1));
	
    //vec3 eye, lookat, up;
    //float fov = camerainfo(eye, lookat, up);
	
    //gl_FragColor = texture2D(tex0, texcoord);
    //gl_FragColor.xy = texcoord.xy;
    //vec3 lt = normalize(vec3(1,1,1));
    vec3 p;
    vec3 n;
    vec3 dir;
    isectinfo(p, n, dir);

    float ndotl = 0.0;
    vec3 nn = faceforward(n, dir, n);

    int depth = 0;
    depth = raydepth(depth);
    float a = dot(n,lightdir);
    a = ( a + 3.0) * 0.25;
    vec3 color = colors * a;
    gl_FragColor = vec4(color, 1.0);
    return;
    
#if 0
    if (depth > 1) {
        // 2nd ray will have gray color(0.25)
        float ra;
        ra = rayattrib(ra);
        gl_FragColor = vec4(ra);
        return;
    } else {
        vec3 pp = p + 0.01 * dir;
        vec4 hitcol;
        float rayattribute = 0.25;
        float dist = trace(pp, dir, hitcol, rayattribute);
        if (dist > -100.0) {
            gl_FragColor = vec4(hitcol);
        } else {
            gl_FragColor = vec4(0.0, 1.0, 0.0, 1.0);
        }
    }
    return;
#endif

#if 1
    vec3 col = matToColor(matID);

    ndotl = dot(nn, normalize(vec3(1.0, 1.0, -1.0)));
    ndotl = max(0.2, ndotl);  // 0.2 = ambient
    //gl_FragColor = vec4(col * ndotl, 1.0);
    //return;
	
    //vec3 nn = n;
    //if (dot(n, dir) > 0.0) {
    //    nn = -nn;
    //}
	
	
    // trace ray
    //int depth;
    //raydepth(depth);

    float ra;
    ra = rayattrib(ra);
	
	vec3 N = normalize(n);
    vec3 pp = p + 0.01 * nn;
	float ao = 0.0;
	float num = 16.0;
	for (float i = 0.0; i < num; i += 1.0) {
		float theta = sqrt(noise(p.xy,i/num*2.0*PI));
		float phi   = 2.0 * PI * noise(p.xy+vec2(0.3,0.2),i/num*2.0*PI+23.5);
		vec3 ref;
		ref.x = cos(phi) * theta;
		ref.y = sin(phi) * theta;
		ref.z = sqrt(1.0 - theta * theta);
		vec3 basis0 = vec3(1,0,0);
		vec3 basis1 = vec3(1,1,1);
		vec3 basis2 = vec3(0,0,1);
		/*basis0 = cross(N,basis1);
		basis2 = cross(basis0,basis1);
		basis1 = cross(basis2,basis0);*/
		orthoBasis(basis0,basis1,basis2, nn);
		vec3 rray;
		rray.x = ref.x * basis0.x + ref.y * basis1.x + ref.z * basis2.x;
		rray.y = ref.x * basis0.y + ref.y * basis1.y + ref.z * basis2.y;
		rray.z = ref.x * basis0.z + ref.y * basis1.z + ref.z * basis2.z;
		
		float dist = trace(pp, rray);//lightdir);
		if (dist < -1000.0) {
			ao += 1.0;
		}
	}
    float ambient = 0.2;
	
	float y = ao / num;
	/*    float y;
	 if (dist < -100.0) {
	 y = (1.0 - max(cos(max(0,dot(-lightdir, normalize(nn))) * 3.141592), ambient)) * 0.5;
	 } else {
	 y = ambient; // shadowed
	 }*/

    ndotl = dot(nn, normalize(vec3(1.0, 1.0, 1.0)));
    ndotl = max(0.2, y * ndotl);  // 0.2 = ambient
	
    //gl_FragColor = vec4(nn, 1);
    gl_FragColor = vec4(ndotl * col, 1.0);
#endif
    gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);
}
