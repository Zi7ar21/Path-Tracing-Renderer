#include <iostream>

//#include <math.h>

// OpenGL Mathematics
#include <glm/glm.hpp>
using namespace glm;

// isinf() and isnan() for Microsoft Visual C++
#ifdef _MSC_VER
#define isinff(a) isinf(a)
#define isnanf(a) _isnanf(a)
#endif

// ##### Parameters #####

// Camera Field of View
#define CAMERA_FOV 0.5f

// Tonemap Exposure Parameter (only if HDR is disabled)
#define EXPOSURE 1.0f

// Portable FloatMap
//#define HDR

// Ray-Marching Tolerance
#define HIT_DIST 0.003f

// Maximum Path-Tracing Bounces
#define MAX_BOUNCES 4u

// Maximum Samples
#define MAX_SAMPLES 32u

// Maximum Ray-Marching Steps
#define MAX_STEPS 512u

// ##### Constants #####

// http://www.mimirgames.com/articles/programming/digits-of-pi-needed-for-floating-point-numbers/
const float pi         = 3.141592653589793f;
const float inverse_pi = 0.318309886183790f;

// ##### Classes/Types #####

class material_properties
{
	public:
		unsigned int material_id;
		vec3         base_color ;
		float        roughness  ;
};

class raycast
{
	public:
		bool  expire;
		bool  hit   ;
		float tMin  ;
		float tMax  ;
		vec3  normal;
		material_properties material_data;
};

class image_buffer
{
	public:
		unsigned int size_x;
		unsigned int size_y;
		unsigned int size;
		vec3 *buffer;
	public:
		void allocate(unsigned int buffer_size_x, unsigned int buffer_size_y)
		{
			size_x = buffer_size_x;
			size_y = buffer_size_y;
			size = buffer_size_x * buffer_size_y;
			buffer = (vec3*)malloc(sizeof(vec3) * size);
		}

		void cleanup()
		{
			free(buffer);
		}
};

// https://nullprogram.com/blog/2018/07/31/
//uint32_t triple32(uint32_t x);

// Random Number Generator Seed
uint32_t ns;

void init_rng(uint32_t seed)
{
	ns = seed;
}

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
float random_float()
{
	// Update RNG
	ns = triple32(ns);

	return float(ns)/float(0xFFFFFFFFu);
}

// 2-Component Uniform Random Vector
vec2 rand2()
{
	vec2 vector;
	vector.x = random_float();
	vector.y = random_float();
	return vector;
}

// 3-Component Uniform Random Vector
vec3 rand3()
{
	vec3 vector;
	vector.x = random_float();
	vector.y = random_float();
	vector.z = random_float();
	return vector;
}

// 4-Component Uniform Random Vector
vec4 rand4()
{
	vec4 vector;
	vector.x = random_float();
	vector.y = random_float();
	vector.z = random_float();
	vector.w = random_float();
	return vector;
}

// See michael0884's usage of PCG Random
// https://www.shadertoy.com/view/wltcRS
// https://www.shadertoy.com/view/WttyWX

vec2 nrand2(vec2 mean, float sigma)
{
	vec2 z = rand2();
	return mean + sigma * glm::sqrt( -2.0f * glm::log(z.x) ) * vec2( glm::cos(2.0f * pi * z.y), glm::sin(2.0f * pi * z.y) );
}

vec3 nrand3(vec3 mean, float sigma)
{
	vec4 z = rand4();
	return mean + sigma * glm::sqrt( -2.0f * glm::log( vec3(z.x, z.x, z.y ) ) ) * vec3( glm::cos(2.0f * pi * z.z), glm::sin(2.0f * pi * z.z), glm::cos(2.0f * pi * z.w) );
}

vec4 nrand4(vec4 mean, float sigma)
{
	vec4 z = rand4();
	return mean + sigma * glm::sqrt( -2.0f * glm::log( vec4(z.x, z.x, z.y, z.y) ) ) * vec4( glm::cos(2.0f * pi * z.z), glm::sin(2.0f * pi * z.z), glm::cos(2.0f * pi * z.w), glm::sin(2.0f * pi * z.w) );
}

// Random Uniform Direction
vec2 udir2()
{
	float z = rand();
	float r = 2.0f * pi * z;
	float s = glm::sin(r), c = glm::cos(r);
	return vec2(s, c);
}

// Random Uniform Direction
vec3 udir3()
{
	vec2 z = rand2();
	vec2 r = vec2( 2.0f * pi * z.x, acosf(2.0f * z.y - 1.0f) );
	vec2 s = glm::sin(r), c = glm::cos(r);
	return vec3(s.y * c.x, s.x * s.y, c.y);
}

// http://blog.hvidtfeldts.net/index.php/2011/09/distance-estimated-3d-fractals-v-the-mandelbulb-different-de-approximations/
float mandelbulb(vec3 p)
{
	const float power = 8;

	p = vec3(p.x, p.z, p.y);

	vec3 z = p;
	float dr = 1;
	float r  = 0;

	for(unsigned int i = 0; i < 128; i++)
	{
		r = length(z);

		if(r > 4.0f)
		{
			break;
		}

		// Convert to Polar Coordinates
		float theta = glm::acos(z.z/r   );
		float phi   = glm::atan(z.y, z.x);
		dr = glm::pow(r, power - 1.0f) * power * dr + 1.0;

		// Scale and Rotate the Point
		float zr = glm::pow(r, power);
		theta = theta * power;
		phi   = phi   * power;

		float sin_theta = glm::sin(theta);
		float cos_theta = glm::cos(theta);
		float sin_phi = glm::sin(phi);
		float cos_phi = glm::cos(phi);

		// Convert back to Cartesian Coordinates
		z = zr * vec3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);

		z += p;
	}

	return 0.5f * glm::log(r) * r / dr;
}

/*
float DE(vec3 p, unsigned int *id)
{
	*id = 1;

	float l = length(p);

	if(l > 2.3f)
	{
		return l - 2.2f;
	}

	return mandelbulb(p);
}
*/

float DE(vec3 p, unsigned int *id)
{
	float DE0 = p.y + 1.0f; // Plane
	//float DE0 = length( p - vec3(-0, -1, -0) ) - 1.0f;
	float DE1 = length(p) - 1.0f; // Sphere

	float minDE = glm::min(DE0, DE1);

	*id = 0; // Null
	*id = (minDE == DE0) ? 1 : *id; // Plane
	*id = (minDE == DE1) ? 2 : *id; // Sphere

	return minDE;
}

void update_material(raycast *raycast_data)
{
	if(raycast_data->material_data.material_id == 1)
	{
		raycast_data->material_data.base_color = vec3(0.8f);
		raycast_data->material_data.roughness  = 0.3f;
	}
	else if(raycast_data->material_data.material_id == 2)
	{
		raycast_data->material_data.base_color = vec3(0.8f);
		raycast_data->material_data.roughness  = 0.3f;
	}
	else
	{
		raycast_data->material_data.base_color = vec3(1.0f, 0.0f, 1.0f);
		raycast_data->material_data.roughness  = 0.3f;
	}
}

vec3 ortho(vec3 v)
{
	return normalize( glm::abs(v.x) > glm::abs(v.z) ? vec3(-v.y, v.x, 0.0f) : vec3(0.0f, -v.z, v.y) );
}

// Weighted Cosine Sampling
vec3 cosine_weighted_sample()
{
	vec2 z = rand2();
	vec2 r = vec2( 2.0f * pi * z.x, glm::sqrt(z.y) );
	return vec3( r.y * vec2( glm::cos(r.x), glm::sin(r.x) ), glm::sqrt(1.0f - r.y * r.y) );
}

// Lambertian Probability Distribution Function
float PDF(vec3 wi, vec3 wo, material_properties material_data)
{
	return glm::max(wi.z, 0.001f) * inverse_pi;
}

// Lambertian Bidirectional Reflectance Distribution Function
vec3 BRDF(vec3 wi, vec3 wo, material_properties material_data)
{
	return material_data.base_color * inverse_pi;
}

vec3 sky_radiance(vec3 dir)
{
	const vec3 sun_dir = normalize( vec3( 1.0f,  1.0f, -1.0f) );

	float d = dot(dir, sun_dir);

	if(d > 0.8f)
	{
		return vec3(5);
	}

	return vec3(1) * (glm::max( (0.8f * dir.y) + 0.15f, 0.0f ) + 0.05f);
}

// SDF Tetrahedron Numerical Normals
vec3 calculate_normal(vec3 p)
{
	/*
	unsigned int null_id;

	return normalize(
	vec3(-1, -1, -1) * DE(p + vec3(-1, -1, -1) * HIT_DIST, &null_id) +
	vec3(-1,  1,  1) * DE(p + vec3(-1,  1,  1) * HIT_DIST, &null_id) +
	vec3( 1, -1,  1) * DE(p + vec3( 1, -1,  1) * HIT_DIST, &null_id) +
	vec3( 1,  1, -1) * DE(p + vec3( 1,  1, -1) * HIT_DIST, &null_id)
	);
	*/

	unsigned int null_id;

	vec3 v0 = vec3(-1.0f, -1.0f, -1.0f) * DE(p + (vec3(-1.0f, -1.0f, -1.0f) * HIT_DIST), &null_id);
	vec3 v1 = vec3(-1.0f,  1.0f,  1.0f) * DE(p + (vec3(-1.0f,  1.0f,  1.0f) * HIT_DIST), &null_id);
	vec3 v2 = vec3( 1.0f, -1.0f,  1.0f) * DE(p + (vec3( 1.0f, -1.0f,  1.0f) * HIT_DIST), &null_id);
	vec3 v3 = vec3( 1.0f,  1.0f, -1.0f) * DE(p + (vec3( 1.0f,  1.0f, -1.0f) * HIT_DIST), &null_id);

	return normalize(v0 + v1 + v2 + v3);
}

bool inside_bounding_box(vec3 p, vec3 minimum, vec3 maximum)
{
	return
	p.x > minimum.x && p.x < maximum.x &&
	p.y > minimum.y && p.y < maximum.y &&
	p.z > minimum.z && p.z < maximum.z;
}

/*
// Axis-Aligned Bounding Box Intersection
bool intersect_aabb(vec3 minimum, vec3 maximum, vec3 origin, vec3 direction, float *tmin, float *tmax)
{
	vec3 inv_direction = 1.0f / direction;
	vec3 tmin_vec = (minimum - origin) * inv_direction;
	vec3 tmax_vec = (maximum - origin) * inv_direction;

	*tmin = max(max(tmin_vec.x, tmin_vec.y), tmin_vec.z);
	*tmax = min(min(tmax_vec.x, tmax_vec.y), tmax_vec.z);

	return *tmin <= *tmax;
}
*/

raycast process_hit(raycast raycast_data, raycast raycast_data_)
{
	if(raycast_data_.hit)
	{
		if(raycast_data.hit)
		{
			raycast output_data;

			output_data.expire = false;
			output_data.hit = true;

			if(raycast_data.tMin < raycast_data_.tMin)
			{
				output_data.tMin          = raycast_data.tMin;
				output_data.normal        = raycast_data.normal;
				output_data.material_data = raycast_data.material_data;
			}
			else
			{
				output_data.tMin          = raycast_data_.tMin;
				output_data.normal        = raycast_data_.normal;
				output_data.material_data = raycast_data_.material_data;
			}

			if(raycast_data.tMax > raycast_data_.tMax)
			{
				output_data.tMax = raycast_data.tMax;
			}
			else
			{
				output_data.tMax = raycast_data_.tMax;
			}

			return output_data;
		}

		return raycast_data_;
	}

	return raycast_data;
}

raycast trace(vec3 ro, vec3 rd)
{
	raycast raycast_data;

	raycast_data.expire = true   ;
	raycast_data.hit    = false  ;
	raycast_data.tMin   = -1.0f  ;
	raycast_data.tMax   = -1.0f  ;
	raycast_data.normal = vec3(0);
	raycast_data.material_data.material_id = 0;

	float t = 0;
	unsigned int id = 0;

	for(unsigned int i = 0; i < MAX_STEPS; i++)
	{
		vec3 rayPos = ro + rd * t;

		//if(glm::abs(rayPos.x) > 2.3f || glm::abs(rayPos.y) > 2.3f || rayPos.z > 4.01f || rayPos.z < -2.3f)
		if(glm::abs(rayPos.x) > 8.0f || glm::abs(rayPos.y) > 8.0f || glm::abs(rayPos.z) > 8.0f)
		{
			raycast_data.expire = false;

			break;
		}

		float td = DE(rayPos, &id);

		if(td < HIT_DIST)
		{
			raycast_data.expire = false;
			raycast_data.hit    = true ;
			raycast_data.tMin   = td   ;
			raycast_data.normal = calculate_normal(rayPos);
			raycast_data.material_data.material_id = id;

			break;
		}

		t += td;
	}

	return raycast_data;
}

vec3 radiance(vec3 ro, vec3 rd)
{
	vec3 rayPos = ro;
	vec3 rayDir = rd;

	vec3 attenuation = vec3(1);

	for(unsigned int bounces = 0; bounces < MAX_BOUNCES; bounces++)
	{
		raycast raycast_data = trace(rayPos, rayDir);

		if(raycast_data.expire)
		{
			break;
		}

		if(!raycast_data.hit)
		{
			return sky_radiance(rayDir) * attenuation;
		}

		float tMin = raycast_data.tMin;
		vec3  n    = raycast_data.normal;
		material_properties m = raycast_data.material_data;

		update_material(&raycast_data);

		vec3 t = ortho(   n);
		vec3 b = cross(t, n);
		mat3 surf2world = mat3(t, b, n);
		mat3 world2surf = transpose(surf2world);

		vec3 wi = cosine_weighted_sample();
		vec3 wo = world2surf * -rayDir;

		// loicvdb's magic wisdom go brrr
		attenuation *= ( BRDF(wi, wo, m) / PDF(wi, wo, m) ) * glm::max(wi.z, 0.0f);

		rayPos = rayPos + ( rayDir * (tMin - HIT_DIST) );
		//rayPos += (rayDir * tMin) + (HIT_DIST * raycast_data.normal);

		rayDir = surf2world * wi;

		rayDir = normalize(rayDir);
	}

	return vec3(-1);
}

// Blackman-Harris Pixel Filter
vec2 pixel_filter(vec2 pixel_coord)
{
	// https://en.wikipedia.org/wiki/Window_function#Blackman–Harris_window
	// w[n] = a0 - a1 * cos(2 * pi * n / N) + a2 * cos(4 * pi * n / N) - a3 * cos(6 * pi * n / N)
	// a0 = 0.35875; a1 = 0.48829; a2 = 0.14128; a3 = 0.01168;

	const float a0 = 0.35875f;
	const float a1 = 0.48829f;
	const float a2 = 0.14128f;
	const float a3 = 0.01168f;

	//float n = 0.5f * random_float() + 0.5f;
	float n = random_float();

	float w = a0 - a1 * glm::cos(2.0f * pi * n) + a2 * glm::cos(4.0f * pi * n) - a3 * glm::cos(6.0f * pi * n);

	return pixel_coord + (2.0f * udir2() * w);
}

#include <fstream>
#include <sstream>

void write_render(image_buffer render_buffer)
{
	std::ofstream image_file;

	image_file.open("render.ppm");

	image_file << "P6\n" << render_buffer.size_x << " " << render_buffer.size_y << "\n255\n";

	for(unsigned int i = 0; i < render_buffer.size; i++)
	{
		// Quantization
		unsigned char channel_r = (unsigned char)std::min(std::max( (unsigned int)(255.0 * render_buffer.buffer[i].r), (unsigned int)0 ), (unsigned int)255);
		unsigned char channel_g = (unsigned char)std::min(std::max( (unsigned int)(255.0 * render_buffer.buffer[i].g), (unsigned int)0 ), (unsigned int)255);
		unsigned char channel_b = (unsigned char)std::min(std::max( (unsigned int)(255.0 * render_buffer.buffer[i].b), (unsigned int)0 ), (unsigned int)255);

		image_file << channel_r << channel_g << channel_b;
	}

	image_file.close();
}

void write_render_HDR(image_buffer render_buffer)
{
	std::ofstream image_file;

	image_file.open("render.pfm");

	image_file << "PF\n" << render_buffer.size_x << " " << render_buffer.size_y << "\n-1.0\n";

	for(unsigned int i = 0; i < render_buffer.size; i++)
	{
		// https://stackoverflow.com/questions/30923685/writing-floats-to-a-binary-file-in-c-equivalent-of-javas-dataoutputstream-w
		image_file.write( reinterpret_cast<const char*>(&render_buffer.buffer[i].r), sizeof(float) );
		image_file.write( reinterpret_cast<const char*>(&render_buffer.buffer[i].g), sizeof(float) );
		image_file.write( reinterpret_cast<const char*>(&render_buffer.buffer[i].b), sizeof(float) );
	}

	image_file.close();
}

void write_frame(image_buffer render_buffer, unsigned int frame_number)
{
	std::stringstream filename; filename.fill('0'); filename.width(3); filename << std::to_string(frame_number);

	std::string file_name = "render/frame" + filename.str() + ".ppm";

	std::ofstream image_file;

	image_file.open(file_name);

	image_file << "P6\n" << render_buffer.size_x << " " << render_buffer.size_y << "\n255\n";

	for(unsigned int i = 0; i < render_buffer.size; i++)
	{
		// Quantization
		unsigned char channel_r = (unsigned char)std::min(std::max( (unsigned int)(255.0 * render_buffer.buffer[i].r), (unsigned int)0 ), (unsigned int)255);
		unsigned char channel_g = (unsigned char)std::min(std::max( (unsigned int)(255.0 * render_buffer.buffer[i].g), (unsigned int)0 ), (unsigned int)255);
		unsigned char channel_b = (unsigned char)std::min(std::max( (unsigned int)(255.0 * render_buffer.buffer[i].b), (unsigned int)0 ), (unsigned int)255);

		image_file << channel_r << channel_g << channel_b;
	}

	image_file.close();
}

/*
// https://www.iquilezles.org/www/articles/cputiles/cputiles.htm
void render_image(image_buffer render_buffer)
{
}
*/

void render_image(image_buffer render_buffer)
{
	// Initialize Random Number Generator
	init_rng(1);

	vec2 resolution = vec2(render_buffer.size_x, render_buffer.size_y);

	for(unsigned int pixel_coord_x = 0; pixel_coord_x < render_buffer.size_x; pixel_coord_x++) {
	for(unsigned int pixel_coord_y = 0; pixel_coord_y < render_buffer.size_y; pixel_coord_y++) {
		vec3 color = vec3(0);
		unsigned int samples = 0;

		for(unsigned int i = 0; i < MAX_SAMPLES; i++)
		{
			//init_rng(i);

			vec2 uv = 2.0f * ( ( pixel_filter( vec2(pixel_coord_x, pixel_coord_y) ) - ( 0.5f * vec2(resolution) ) ) / glm::max(resolution.x, resolution.y) );

			vec3 ro = vec3(0.0f, 0.0f, 4.0f);
			vec3 rd = normalize( vec3(CAMERA_FOV * uv, -1) );

			vec3 c = radiance(ro, rd);

			bool sample_expired = c.r < 0.0f || c.g < 0.0f || c.b < 0.0f;
			bool sample_invalid = isinff(c.r) || isinff(c.g) || isinff(c.b) || isnanf(c.r) || isnanf(c.g) || isnanf(c.b);

			// Check if the sample was discarded and accumulate
			if(!sample_expired && !sample_invalid)
			{
				color += c;

				samples++;
			}
		}

		color = samples != 0 ? color / float(samples) : color;

		#ifdef HDR
		render_buffer.buffer[pixel_coord_x + (pixel_coord_y * render_buffer.size_x)] = color;
		#else
		color = clamp(1.0f - glm::exp(-glm::max(color, 0.0f) * EXPOSURE), 0.0f, 1.0f);
		render_buffer.buffer[pixel_coord_x + ( ( (render_buffer.size_y - 1) - pixel_coord_y ) * render_buffer.size_x )] = color;
		#endif
	}
	}
}

// ##### Main #####
int main()
{
	const unsigned int RENDER_SIZE_X = 256;
	const unsigned int RENDER_SIZE_Y = 192;

	std::cout << "Initializing..." << std::endl;

	image_buffer render_buffer;

	render_buffer.allocate(RENDER_SIZE_X, RENDER_SIZE_Y);

	std::cout << "Starting Render..." << std::endl;

	render_image(render_buffer);

	std::cout << "Writing Render to Disk..." << std::endl;

	#ifdef HDR
	write_render_HDR(render_buffer);
	#else
	write_render(render_buffer);
	#endif

	std::cout << "Cleaning Up..." << std::endl;

	render_buffer.cleanup();

	std::cout << "Done!" << std::endl;

	return EXIT_SUCCESS;
}