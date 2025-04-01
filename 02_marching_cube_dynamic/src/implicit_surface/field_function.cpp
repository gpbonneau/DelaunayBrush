#include "field_function.hpp"
#include <algorithm> // For std::max_element
#include <vector>    // For std::vector

using namespace cgp;

// Parameterization of a sphere (using exact distance field)
static float sphere_distance(vec3 const& p, vec3 const& center, float radius)
{
    return norm(p - center) - radius; // Exact signed distance field
}

float field_function_structure::operator()(cgp::vec3 const& p) const
{
    if (spheres.empty()) {
        return 1.0f; // Return background field if no spheres
    }

    // Compute distance to each sphere
    std::vector<float> distances;
    distances.reserve(spheres.size());
    
    for (const auto& sphere : spheres) {
        distances.push_back(sphere_distance(p, sphere.center, sphere.radius));
    }

    // Boolean union operation (minimum of distances)
    float value = *std::min_element(distances.begin(), distances.end());

    // Optional: Add noise to the final result
    if (noise_magnitude > 0) {
        vec3 const offset = vec3{ noise_offset + 1000, 1000, 1000 };
        vec3 const p_noise = noise_scale * p + offset;
        value += noise_magnitude * noise_perlin(p_noise, noise_octave, noise_persistance);
    }

    return value;
}