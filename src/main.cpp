#include <iostream>

#include <common.hpp>

#include <fstream>
#include <sstream>

#define HIT_DIST 0.003f

#define CAMERA_FOV 1.0f

float dot_p(vec2 vector) { return dot(vector, vector); }
float dot_p(vec3 vector) { return dot(vector, vector); }
float dot_p(vec4 vector) { return dot(vector, vector); }

float map(vec3 p)
{
	float d0 = dot_p( p - vec3(-1.0f, -1.0f, -1.0f) );
	float d1 = dot_p( p - vec3( 1.0f, -1.0f, -1.0f) );
	float d2 = dot_p( p - vec3(-1.0f,  1.0f, -1.0f) );
	float d3 = dot_p( p - vec3( 1.0f,  1.0f, -1.0f) );
	float d4 = dot_p( p - vec3(-1.0f, -1.0f,  1.0f) );
	float d5 = dot_p( p - vec3( 1.0f, -1.0f,  1.0f) );
	float d6 = dot_p( p - vec3(-1.0f,  1.0f,  1.0f) );
	float d7 = dot_p( p - vec3( 1.0f,  1.0f,  1.0f) );

	float min_d = glm::min( glm::min( glm::min(d0, d1), glm::min(d2, d3) ), glm::min( glm::min(d4, d5), glm::min(d6, d7) ) );

	min_d = glm::sqrt(min_d) - 0.5f;

	min_d = glm::min(min_d, p.y + 1.0f);

	return min_d;
}

vec3 calculate_normal(vec3 p)
{
	return normalize(
	vec3(-1, -1, -1) * map(p + vec3(-1, -1, -1) * HIT_DIST) +
	vec3(-1,  1,  1) * map(p + vec3(-1,  1,  1) * HIT_DIST) +
	vec3( 1, -1,  1) * map(p + vec3( 1, -1,  1) * HIT_DIST) +
	vec3( 1,  1, -1) * map(p + vec3( 1,  1, -1) * HIT_DIST)
	);
}

float trace(vec3 ro, vec3 rd)
{
	float t = 0.0f;

	for(unsigned int i = 0u; i < 128u; i++)
	{
		vec3 p = ro + rd * t;

		if(glm::abs(p.x) > 4.0f || glm::abs(p.y) > 4.0f || glm::abs(p.z) > 4.0f)
		{
			break;
		}

		float d = map(p);

		if(d < HIT_DIST)
		{
			return t;
		}

		t += d;
	}

	return -1.0f;
}

vec3 sky_radiance(vec3 dir)
{
	if(dir.y > 0.8f)
	{
		return vec3(1);
	}

	return vec3(0);
}

vec3 radiance(vec3 ro, vec3 rd, uint32_t *rng_state)
{
	vec3 ray_pos = ro;
	vec3 ray_dir = rd;

	vec3 att = vec3(1);

	for(unsigned int bounces = 0u; bounces < 4u; bounces++)
	{
		float t = trace(ray_pos, ray_dir);

		if(t < 0.0f)
		{
			return att * sky_radiance(ray_dir);
		}

		vec3 n = calculate_normal(ro + rd * t);

		ray_pos += ray_dir * (t - HIT_DIST);
		ray_dir = reflect( ray_dir, normalize( nrand3(n, 1.0f, rng_state) ) );

		att *= vec3(0.8f, 0.8f, 0.8f);
	}

	return vec3(-1);
}

void render_pixel(vec3 *pixel, uvec2 resolution, uvec2 pixel_coord)
{
	uint32_t ns = pixel_coord.x + pixel_coord.y * resolution.x + 1u;

	vec3 acc = vec3(0);
	unsigned int s = 0u;

	for(unsigned int samples = 0u; samples < 128u; samples++)
	{
		vec2 uv = 2.0f * ( nrand2(vec2(pixel_coord), 0.5f, &ns) - 0.5f * vec2(resolution) ) / float( glm::max(resolution.x, resolution.y) );

		vec3 ro = vec3(0.0f, 0.0f, 2.0f);
		vec3 rd = normalize( vec3(CAMERA_FOV * uv, -1.0f) );

		vec3 c = radiance(ro, rd, &ns);

		if(c.r >= 0.0f && c.g >= 0.0f && c.b >= 0.0f)
		{
			acc += c;
			s++;
		}
	}

	acc = s != 0u ? acc / float(s) : acc;

	*pixel = acc;
}

void render_image(vec3 *image, uvec2 resolution)
{
	/*
	uvec2 pixel_coord;
	for(pixel_coord.x = 0u; pixel_coord.x < resolution.x; pixel_coord.x++) {
	for(pixel_coord.y = 0u; pixel_coord.y < resolution.y; pixel_coord.y++) {
		render_pixel(&image[pixel_coord.x + pixel_coord.y * resolution.x], resolution, pixel_coord);
	}
	}
	*/
	unsigned int TILE_SIZE = 32;

	// prep data
	unsigned int numxtiles = resolution.x / TILE_SIZE;
	unsigned int numytiles = resolution.y / TILE_SIZE;
	unsigned int numtiles  = numxtiles * numytiles;

	// render tiles
	#pragma omp parallel for
	for(unsigned int tile = 0u; tile < numtiles; tile++)
	{
		// tile offset
		unsigned int ia = TILE_SIZE * (tile % numxtiles);
		unsigned int ja = TILE_SIZE * (tile / numxtiles);

		// for every pixel in this tile, compute color
		for(unsigned int j = 0u; j < TILE_SIZE; j++) {
		for(unsigned int i = 0u; i < TILE_SIZE; i++) {
			render_pixel( &image[resolution.x*(ja+j)+(ia+i)], resolution, uvec2(ia+i, ja+j) );
		}
		}
	}

}

void write_render_HDR(vec3 *image, uvec2 resolution)
{
	std::ofstream image_file;

	image_file.open("render.pfm");

	image_file << "PF\n" << resolution.x << " " << resolution.y << "\n-1.0\n";

	for(unsigned int i = 0u; i < resolution.x * resolution.y; i++)
	{
		// https://stackoverflow.com/questions/30923685/writing-floats-to-a-binary-file-in-c-equivalent-of-javas-dataoutputstream-w
		image_file.write( reinterpret_cast<const char*>(&image[i].r), sizeof(float) );
		image_file.write( reinterpret_cast<const char*>(&image[i].g), sizeof(float) );
		image_file.write( reinterpret_cast<const char*>(&image[i].b), sizeof(float) );
	}

	image_file.close();
}

int main()
{
	uvec2 resolution = uvec2(512u, 256u);

	vec3 *image;

	image = (vec3*)malloc(sizeof(vec3) * resolution.x * resolution.y);

	render_image(image, resolution);

	write_render_HDR(image, resolution);

	free(image);

	std::cout << "Done!" << std::endl;

	return EXIT_SUCCESS;
}