#pragma once

#include "cgp/cgp.hpp"

// Parametric function defined as a sum of blobs-like primitives
//  f(p) = sa exp(-||p-pa||^2) + sb exp(-||p-pb||^2) + sc exp(-||p-pc||^2) + noise(p)
//   with noise: a Perlin noise
// The operator()(vec3 p) allows to query a value of the function at arbitrary point in space

struct Sphere {
    cgp::vec3 center;
    double radius;
};

struct field_function_structure
{

	// Query the value of the function at any point p
	float operator()(cgp::vec3 const &p) const;

	// ***************************//
	// Parameters
	// ***************************//

	// The 3 center of the 3 blobs
	cgp::vec3 a = {1.0f, 0.0f, 0.0f};
	cgp::vec3 b = {-25.0f, 20.5f, 0.5f};
	cgp::vec3 c = {0.0F, 1.0f, 0.0f};
	cgp::vec3 d = {0.0F, 0.0f, 1.0f};
	cgp::vec3 e = {2.0F, 1.0f, 0.0f};
	std::vector<Sphere> spheres;

	// The magnitude associated to each primitive
	float sa = 1.0f;
	float sb = 1.0f;
	float sc = 0.0f; // note: The third blob is not visible initially has its magnitude is 0

	// The parameters of the Perlin noise
	float noise_magnitude = 0.0f;	// Magnitude of the noise
	float noise_offset = 0.0f;		// An offset in the parametric domain (get a different value of noise with same parameters)
	float noise_scale = 1.0f;		// Scale in the parametric domain
	int noise_octave = 5;			// Maximum number of octave
	float noise_persistance = 0.3f; // Persistence of Perlin noise
};
