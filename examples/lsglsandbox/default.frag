#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
uniform vec2 resolution;

const int ITER = 300;
const float CHOP = 10.;
const vec2 CENTER = vec2(0., -0.65);
const float bpm = 130.;
const float loop = 60. / bpm * 16.;

float diverge(vec2 c)
{
	vec2 z = vec2(0., 0.);

	for (int i = 0; i < ITER; ++i) {
		z = vec2(z.x * z.x - z.y * z.y, 2. * z.x * z.y) + c;
		if (length(z) > 2.) {
			return float(i);
		}
	}
	return -1.;
}

vec3 color(float count, vec2 pos)
{
	float k = 10.;
	if (count < 0.) {
		if (mod(time, loop) <= loop / 2.) {
			return vec3(0., 0., 0.);
		} else {
			if (mod(floor(k * pos.x) + floor(k * pos.y), 2.) > 0.5) {
				return vec3(1., 0., 0.);
			} else {
				return vec3(1., 1., 1.);
			}
		}
	}
	return vec3(1., mod(count, CHOP) / CHOP, 0.);
}

float delta()
{
	float spb = 60. / bpm;
	float rate = 0.25;
	float beat_pos = mod(time, spb) / spb;
	if (0. <= beat_pos && beat_pos <= rate) {
		return 1. + (1. - beat_pos / rate) * 0.2;
	} else {
		return 1. + (beat_pos - rate) / (1. - rate) * 0.2;
	}
}

float flash()
{
	float spb = 60. / bpm;
	return mod(time, spb) / spb * 0.5 + 0.5;
}

vec2 rotate(vec2 src)
{
	float k = 1.;
	src -= CENTER;
	vec2 result = vec2(src.x * sin(k * time) + src.y * cos(k * time), src.x * cos(k * time) - src.y * sin(k * time));
	result += CENTER;
	return result;
}

bool rect(float sx, float sy, float w, float h, vec2 p)
{
	return sx <= p.x && p.x <= sx + w && sy <= p.y && p.y <= sy + h;
}

bool f(vec2 p)
{
	return rect(0., -0.15, 0.02, 0.15, p) || rect(0., 0., 0.12, 0.02, p) || rect(0., -0.06, 0.10, 0.02, p);
}

bool t(vec2 p)
{
	return rect(0.01, 0., 0.13, 0.02, p) || rect(0.065, -0.15, 0.02, 0.15, p);
}

bool d(vec2 p)
{
	float l = length(p - vec2(0, -0.065));
	return (p.x >= 0. && 0.065 <= l && l <= 0.085) || rect(0., -0.15, 0.02, 0.17, p);
}

bool str(vec2 p)
{
	return t(p-vec2(-0.30, 0.05)) || d(p-vec2(-0.12, 0.05)) || f(p-vec2(0., 0.05));
}

void main(void)
{
	float scale = 1. * (pow(1. - abs((mod(time, loop) / loop - 0.5) * 2.), 1.0) * delta());
	vec2 pos = vec2(-1., 1.) * ((gl_FragCoord.xy / resolution.xy) * 2. - 1.) * resolution / resolution.y * scale + CENTER;
	vec2 rotated = rotate(pos);

	vec3 base_color = color(diverge(rotated), rotated);
	float flash_rate = flash();
	if (mod(time, loop) <= loop / 2. && str(rotated)) {
		gl_FragColor = vec4(1., 1., 1., 1.);
	} else {
		gl_FragColor = vec4(base_color * flash_rate + vec3(1.) * (1. - flash_rate), 1.);
	}
}

