#include "implicit_surface.hpp"
#include <CGAL/Marching_cubes_3.h>
#include <CGAL/Surface_mesh.h>

using namespace cgp;
using Mesh = CGAL::Surface_mesh<CGAL::Point_3<CGAL::Exact_predicates_inexact_constructions_kernel>>;

// Grid resolution and bounds
const int GRID_RES = 50;
const vec3 BOUNDS_MIN(-2, -2, -2);
const vec3 BOUNDS_MAX(2, 2, 2);





static void update_normals(std::vector<vec3>& normals, int number_of_vertex, grid_3D<vec3> const& gradient, std::vector<marching_cube_relative_coordinates> const& relative_coords)
{
	// Compute the normal using linear interpolation of the gradients
	for (int k = 0; k < number_of_vertex; ++k)
	{
		size_t const idx0 = relative_coords[k].k0;
		size_t const idx1 = relative_coords[k].k1;

		vec3 const& n0 = gradient.at_unsafe(idx0);
		vec3 const& n1 = gradient.at_unsafe(idx1);

		float const alpha = relative_coords[k].alpha;
		normals[k] = -normalize((1 - alpha) * n0 + alpha * n1, { 1,0,0 });
	}
}

void implicit_surface_structure::update_marching_cube(float isovalue)
{
	// Variable shortcut
	std::vector<vec3>& position = data_param.position;
	std::vector<vec3>& normal = data_param.normal;
	size_t& number_of_vertex = data_param.number_of_vertex;
	spatial_domain_grid_3D const& domain = field_param.domain;
	std::vector<cgp::marching_cube_relative_coordinates> relative_coord = data_param.relative;
	grid_3D<float> const& field = field_param.field;
	grid_3D<vec3> const& gradient = field_param.gradient;

	// Store the size of the previous position buffer
	size_t const previous_size = position.size();

	// Compute the Marching Cube
	number_of_vertex = marching_cube(position, field.data.data, domain, isovalue, &relative_coord);

	// Resize the vector of normals if needed
	if (normal.size() < position.size())
		normal.resize(position.size());

	update_normals(normal, number_of_vertex, gradient, relative_coord);

	// Update the display of the mesh
	if (position.size() > previous_size) {
		// If there is more position than allocated - perform a full clear and reallocation from scratch
		drawable_param.shape.clear();
		drawable_param.shape.initialize_data_on_gpu(position, normal);
	}
	else {
		// Otherwise simply update the new relevant values re-using the allocated buffers
		drawable_param.shape.vbo_position.update(position, number_of_vertex);
		drawable_param.shape.vbo_normal.update(normal, number_of_vertex);
		drawable_param.shape.vertex_number = number_of_vertex;
	}

}



// Generate a scalar field for Marching Cubes
std::vector<float> generate_scalar_field(field_function_structure& field, int resolution, vec3 min_bound, vec3 max_bound) {
    std::vector<float> field_values(resolution * resolution * resolution);
    
    for (int i = 0; i < resolution; ++i) {
        for (int j = 0; j < resolution; ++j) {
            for (int k = 0; k < resolution; ++k) {
                vec3 p = min_bound + vec3(
                    (max_bound.x - min_bound.x) * i / (resolution - 1),
                    (max_bound.y - min_bound.y) * j / (resolution - 1),
                    (max_bound.z - min_bound.z) * k / (resolution - 1)
                );
                field_values[i * resolution * resolution + j * resolution + k] = field(p);
            }
        }
    }
    return field_values;
}

// Extract mesh using Marching Cubes
Mesh extract_mesh(field_function_structure& field) {
    auto field_values = generate_scalar_field(field, GRID_RES, BOUNDS_MIN, BOUNDS_MAX);
    Mesh mesh;
    
    CGAL::Marching_cubes_3<decltype(field_values.begin())> mc(
        field_values.begin(), field_values.end(), GRID_RES, GRID_RES, GRID_RES,
        BOUNDS_MIN.x, BOUNDS_MIN.y, BOUNDS_MIN.z,
        BOUNDS_MAX.x, BOUNDS_MAX.y, BOUNDS_MAX.z,
        &mesh
    );
    mc.run();
    return mesh;
}

// Save the extracted mesh as an OBJ file
void save_mesh_as_obj(const Mesh& mesh, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    std::unordered_map<Mesh::Vertex_index, int> vertex_index;
    int index = 1;
    for (auto v : mesh.vertices()) {
        auto p = mesh.point(v);
        out << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
        vertex_index[v] = index++;
    }
    
    for (auto f : mesh.faces()) {
        out << "f";
        for (auto v : CGAL::vertices_around_face(mesh.halfedge(f), mesh)) {
            out << " " << vertex_index[v];
        }
        out << "\n";
    }
    out.close();
    std::cout << "Mesh saved to " << filename << std::endl;
}


void implicit_surface_structure::update_field(field_function_structure const& field_function, float isovalue)
{
	// Variable shortcut
	// grid_3D<float>& field = field_param.field;
	// grid_3D<vec3>& gradient = field_param.gradient;
	// spatial_domain_grid_3D& domain = field_param.domain;
	// std::cout<<"update_field"<<std::endl;
	// // Compute the scalar field
	// field = compute_discrete_scalar_field(domain, field_function);
	// std::cout<<"computed discreete"<<std::endl;
	// // Compute the gradient of the scalar field
	// gradient = compute_gradient(field);
	// std::cout<<"computed gradient"<<std::endl;
	// // Recompute the marching cube
	// update_marching_cube(isovalue);
	// std::cout<<"updated_marching_cube"<<std::endl;
	// // Reset the domain visualization (lightweight - can be cleared at each call)
	// drawable_param.domain_box.clear();
	// drawable_param.domain_box.initialize_data_on_gpu(domain.export_segments_for_drawable_border());
	Mesh mesh = extract_mesh(field);
    save_mesh_as_obj(mesh, "output.obj");
}

void implicit_surface_structure::set_domain(int samples, cgp::vec3 const& length)
{
	field_param.domain = spatial_domain_grid_3D::from_center_length({ 0,0,0 }, length, samples * int3{ 1,1,1 });
}



void implicit_surface_structure::gui_update(gui_parameters& gui, field_function_structure& field_function)
{
	bool is_update_marching_cube = false;
	bool is_update_field = false;
	bool is_save_obj = false;

	display_gui_implicit_surface(is_update_field, is_update_marching_cube, is_save_obj, gui, field_function);

	if (is_update_marching_cube) 
		update_marching_cube(gui.isovalue);
	if (is_update_field) {
		set_domain(gui.domain.samples, gui.domain.length);
		update_field(field_function, gui.isovalue);
	}

	if (is_save_obj) {
		data_param.position.resize(data_param.number_of_vertex);
		data_param.normal.resize(data_param.number_of_vertex);
		save_file_obj("mesh.obj", data_param.position, data_param.normal);
	}
		
}



grid_3D<float> compute_discrete_scalar_field(spatial_domain_grid_3D const& domain, field_function_structure const& func)
{
	grid_3D<float> field;
	field.resize(domain.samples);

	// Fill the discrete field values
	for (int kz = 0; kz < domain.samples.z; kz++) {
		for (int ky = 0; ky < domain.samples.y; ky++) {
			for (int kx = 0; kx < domain.samples.x; kx++) {

				vec3 const p = domain.position({ kx, ky, kz });
				field(kx, ky, kz) = func(p);

			}
		}
	}

	return field;
}



grid_3D<vec3> compute_gradient(grid_3D<float> const& field)
{
	grid_3D<vec3> gradient;
	gradient.resize(field.dimension);

	int const Nx = field.dimension.x;
	int const Ny = field.dimension.y;
	int const Nz = field.dimension.z;

	// Use simple finite differences
	//  g(k) = g(k+1)-g(k) // for k<N-1
	//  otherwise g(k) = g(k)-g(k-1)

	for (int kz = 0; kz < Nz; ++kz) {
		for (int ky = 0; ky < Ny; ++ky) {
			for (int kx = 0; kx < Nx; ++kx) {

				vec3& g = gradient.at_unsafe(kx, ky, kz);
				float const f = field.at_unsafe(kx, ky, kz);

				g.x = kx != Nx - 1 ? field.at_unsafe(kx + 1, ky, kz) - f : f - field.at_unsafe(kx - 1, ky, kz);
				g.y = ky != Ny - 1 ? field.at_unsafe(kx, ky + 1, kz) - f : f - field.at_unsafe(kx, ky - 1, kz);
				g.z = kz != Nz - 1 ? field.at_unsafe(kx, ky, kz + 1) - f : f - field.at_unsafe(kx, ky, kz - 1);

			}
		}
	}

	return gradient;
}
