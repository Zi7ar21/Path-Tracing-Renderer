#include <iostream>

#include <common.hpp>

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

/*
// https://www.iquilezles.org/www/articles/cputiles/cputiles.htm
void render_image(image_buffer render_buffer)
{
	const unsigned int TILE_SIZE = 16u;

	// prep data
	const unsigned int numxtiles = resolution.x / TILE_SIZE;
	const unsigned int numytiles = resolution.y / TILE_SIZE;
	const unsigned int numtiles  = numxtiles * numytiles;

	// render tiles
	for(unsigned int tile = 0u; i < numtiles; i++)
	{
		// tile offset
		const unsigned int ia = TILE_SIZE * (tile % numxtiles);
		const unsigned int ja = TILE_SIZE * (tile / numxtiles);

		// for every pixel in this tile, compute color
		for(unsigned int j = 0u; j < TILE_SIZE; j++) {
		for(unsigned int i = 0u; i < TILE_SIZE; i++) {
			image[xres * (ja + j) + (ia + i)] = calcPixelColor(ivec2(ia + i, ja + j), resolution);
		}
		}
	}
}
*/

void render_image(image_buffer render_buffer)
{
	// Initialize Random Number Generator
	init_rng(1u);

	vec2 resolution = vec2(render_buffer.size_x, render_buffer.size_y);

	for(unsigned int pixel_coord_x = 0u; pixel_coord_x < render_buffer.size_x; pixel_coord_x++) {
	for(unsigned int pixel_coord_y = 0u; pixel_coord_y < render_buffer.size_y; pixel_coord_y++) {
		vec3 color = vec3(0.0f);
		unsigned int samples = 0u;

		for(unsigned int i = 0u; i < MAX_SAMPLES; i++)
		{
			//init_rng(i);

			vec2 uv = 2.0f * ( ( pixel_filter( vec2(pixel_coord_x, pixel_coord_y) ) - ( 0.5f * vec2(resolution) ) ) / glm::max(resolution.x, resolution.y) );

			vec3 ro = vec3(0.0f, 0.0f, 4.0f);
			vec3 rd = normalize( vec3(CAMERA_FOV * uv, -1.0f) );

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

		color = samples != 0u ? color / float(samples) : color;

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
	const unsigned int RENDER_SIZE_X = 256u;
	const unsigned int RENDER_SIZE_Y = 192u;

	std::cout << "\nThe Order of the Simulation Path-Tracing Renderer" << std::endl;

	// Output Render Parameters
	std::cout << "\nSettings:\n"
	<< "\nRender Size: " << RENDER_SIZE_X << "x" << RENDER_SIZE_Y
	<< "\nHDR Output: " <<
	#ifdef HDR
	"Yes"
	#else
	"No\nExposure: " << EXPOSURE
	#endif
	<< "\nCamera FOV: " << CAMERA_FOV
	<< "\nMaximum Samples: " << MAX_SAMPLES
	<< "\nMaximum Bounces: " << MAX_BOUNCES
	<< "\nMaximum Ray-Marching Steps: " << MAX_STEPS
	<< "\nRay-Marching Tolerance: " << HIT_DIST << "\n" << std::endl;

	std::cout << "Initializing..." << std::endl;

	image_buffer render_buffer;

	render_buffer.allocate(RENDER_SIZE_X, RENDER_SIZE_Y);

	std::cout << "\nStarting Render..." << std::endl;

	render_image(render_buffer);

	std::cout << "\nWriting Render to Disk..." << std::endl;

	#ifdef HDR
	write_render_HDR(render_buffer);
	#else
	write_render(render_buffer);
	#endif

	std::cout << "\nCleaning Up..." << std::endl;

	render_buffer.cleanup();

	std::cout << "\nDone!\n" << std::endl;

	return EXIT_SUCCESS;
}