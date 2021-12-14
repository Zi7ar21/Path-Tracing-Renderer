#include <common.hpp>

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
	return glm::max(wi.z, 0.000001f) * inverse_pi;
}

// Lambertian Bidirectional Reflectance Distribution Function
vec3 BRDF(vec3 wi, vec3 wo, material_properties material_data)
{
	return material_data.base_color * inverse_pi;
}

vec3 radiance(vec3 ro, vec3 rd)
{
	vec3 rayPos = ro;
	vec3 rayDir = rd;

	vec3 attenuation = vec3(1.0f);

	for(unsigned int bounces = 0u; bounces < MAX_BOUNCES; bounces++)
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

		rayPos += rayDir * (tMin - HIT_DIST);
		//rayPos += (rayDir * tMin) + (HIT_DIST * raycast_data.normal);

		rayDir = surf2world * wi;

		rayDir = normalize(rayDir);
	}

	return vec3(-1.0f);
}