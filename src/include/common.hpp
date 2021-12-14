#include <glm/glm.hpp>
using namespace glm;

// ##### Random Number Generator #####

uint32_t triple32(uint32_t x); // https://nullprogram.com/blog/2018/07/31/

float uniform_random(uint32_t *rng_state); // Random Value Between 0.0 and 1.0

vec2 rand2(uint32_t *rng_state); // 2-Component Uniform Random Vector
vec3 rand3(uint32_t *rng_state); // 3-Component Uniform Random Vector
vec4 rand4(uint32_t *rng_state); // 4-Component Uniform Random Vector


vec2 udir2(uint32_t *rng_state); // Uniformly distributed random point on a unit circle
vec3 udir3(uint32_t *rng_state); // Uniformly distributed random point on the surface of a unit sphere

vec2 nrand2(vec2 mean, float sigma, uint32_t *rng_state);
vec3 nrand3(vec3 mean, float sigma, uint32_t *rng_state);
vec4 nrand4(vec4 mean, float sigma, uint32_t *rng_state);