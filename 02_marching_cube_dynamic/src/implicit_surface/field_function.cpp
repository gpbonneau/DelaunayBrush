#include "field_function.hpp"

using namespace cgp;

// Parameterization Gaussian centered at point p0 with radius sigma
static float gaussian(vec3 const& p, vec3 const& p0, float sigma)
{
    float const d = norm(p - p0);
    float const value = std::exp(-(d * d) / (sigma * sigma));
    return value;
}

// Updated field function to use multiple sphere centers and radii
float field_function_structure::operator()(cgp::vec3 const& p) const
{
    float value = 0.0f;
    // std::cout << p << std::endl;
    // Iterate over all sphere centers and radii
    for (size_t i = 0; i < spheres.size(); ++i)
    {
        value += 1 * gaussian(p, spheres[i].center, spheres[i].radius/2);
    }

    // Add noise contribution if enabled
    // if (noise_magnitude > 0) {
    //     vec3 const offset = vec3{ noise_offset + 1000, 1000, 1000 };
    //     vec3 const p_noise = noise_scale * p + offset;
    //     value += noise_magnitude * noise_perlin(p_noise, noise_octave, noise_persistance);
    // }

    return value;
}
