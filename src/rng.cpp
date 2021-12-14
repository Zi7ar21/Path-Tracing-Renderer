#include <common.hpp>

// https://nullprogram.com/blog/2018/07/31/
uint32_t triple32(uint32_t x)
{
	x ^= x >> 17u;
	x *= 0xED5AD4BBu;
	x ^= x >> 11u;
	x *= 0xAC4C1B51u;
	x ^= x >> 15u;
	x *= 0x31848BABu;
	x ^= x >> 14u;
	return x;
}

// Random Value Between 0.0 and 1.0
float uniform_random(uint32_t *rng_state)
{
	*rng_state = triple32(*rng_state);

	return float(*rng_state) / float(0xFFFFFFFFu);
}

// 2-Component Uniform Random Vector
vec2 rand2(uint32_t *rng_state)
{
	vec2 vector;
	vector.x = uniform_random(rng_state);
	vector.y = uniform_random(rng_state);
	return vector;
}

// 3-Component Uniform Random Vector
vec3 rand3(uint32_t *rng_state)
{
	vec3 vector;
	vector.x = uniform_random(rng_state);
	vector.y = uniform_random(rng_state);
	vector.y = uniform_random(rng_state);
	return vector;
}

// 4-Component Uniform Random Vector
vec4 rand4(uint32_t *rng_state)
{
	vec4 vector;
	vector.x = uniform_random(rng_state);
	vector.y = uniform_random(rng_state);
	vector.z = uniform_random(rng_state);
	vector.w = uniform_random(rng_state);
	return vector;
}

// Uniformly distributied random point on a unit circle
vec2 udir2(uint32_t *rng_state)
{
	float z = uniform_random(rng_state);
	float r = 2.0f * M_PI * z;
	float s = glm::sin(r), c = glm::cos(r);
	return vec2(s, c);
}

// Uniformly distributed random point on the surface of a unit sphere
vec3 udir3(uint32_t *rng_state)
{
	vec2 z = rand2(rng_state);
	vec2 r = vec2( 2.0f * M_PI * z.x, glm::acos(2.0f * z.y - 1.0f) );
	vec2 s = glm::sin(r), c = glm::cos(r);
	return vec3(c.x * s.y, s.x * s.y, c.y);
}

// See michael0884's usage of PCG Random
// https://www.shadertoy.com/view/wltcRS
// https://www.shadertoy.com/view/WttyWX

vec2 nrand2(vec2 mean, float sigma, uint32_t *rng_state)
{
	vec2 z = rand2(rng_state);
	vec2 v = vec2(glm::cos(2.0f * M_PI * z.y), glm::sin(2.0f * M_PI * z.y) );
	return mean + sigma * glm::sqrt( -2.0f * glm::log(z.x) ) * v;
}

vec3 nrand3(vec3 mean, float sigma, uint32_t *rng_state)
{
	vec4 z = rand4(rng_state);
	vec3 v = vec3(glm::cos(2.0f * M_PI * z.z), glm::sin(2.0f * M_PI * z.z), glm::cos(2.0f * M_PI * z.w) );
	return mean + sigma * glm::sqrt( -2.0f * glm::log( vec3(z.x, z.x, z.y) ) ) * v;
}

vec4 nrand4(vec4 mean, float sigma, uint32_t *rng_state)
{
	vec4 z = rand4(rng_state);
	vec4 v = vec4(glm::cos(2.0f * M_PI * z.z), glm::sin(2.0f * M_PI * z.z), glm::cos(2.0f * M_PI * z.w), glm::sin(2.0f * M_PI * z.w) );
	return mean + sigma * glm::sqrt( -2.0f * glm::log( vec4(z.x, z.x, z.y, z.y) ) ) * v;
}