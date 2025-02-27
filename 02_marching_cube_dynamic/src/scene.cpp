#include "scene.hpp"
#include "happly.h"


using namespace cgp;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <igl/readPLY.h>
#include <igl/writePLY.h>
#include <igl/copyleft/cgal/convex_hull.h>
// #include <igl/opengl/glfw/Viewer.h>
#include <unordered_set>
#include <unordered_map>
#include <igl/writeOBJ.h>

// Define CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
typedef K::Point_3 Point;
typedef K::Tetrahedron_3 Tetrahedron;
// typedef CGAL::Tetrahedron_3<Kernel> Tetrahedron;



void generate_sphere_mesh(
    double radius, int resolution,
    Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    std::vector<Eigen::RowVector3d> vertices;
    std::vector<Eigen::RowVector3i> faces;

    // Create vertices
    for (int i = 0; i <= resolution; ++i)
    {
        double theta = M_PI * i / resolution; // Angle from top to bottom
        for (int j = 0; j <= resolution; ++j)
        {
            double phi = 2 * M_PI * j / resolution; // Angle around the sphere
            double x = radius * sin(theta) * cos(phi);
            double y = radius * sin(theta) * sin(phi);
            double z = radius * cos(theta);
            vertices.emplace_back(x, y, z);
        }
    }

    // Create faces
    for (int i = 0; i < resolution; ++i)
    {
        for (int j = 0; j < resolution; ++j)
        {
            int v0 = i * (resolution + 1) + j;
            int v1 = v0 + 1;
            int v2 = v0 + resolution + 1;
            int v3 = v2 + 1;

            faces.emplace_back(v2, v1, v0);
            faces.emplace_back(v2, v3, v1);
        }
    }

    // Convert to Eigen matrices
    V.resize(vertices.size(), 3);
    F.resize(faces.size(), 3);

    for (size_t i = 0; i < vertices.size(); ++i)
        V.row(i) = vertices[i];

    for (size_t i = 0; i < faces.size(); ++i)
        F.row(i) = faces[i];
}
bool readPLYWithNormals(const std::string& ply_file, Eigen::MatrixXd& V, Eigen::MatrixXd& N) {
    // Create a happly PLY data object
    happly::PLYData plyIn(ply_file);

    // Check if the file contains vertex data
    // if (!plyIn.hasElement("vertex") || !plyIn.hasProperty("x") || !plyIn.hasProperty("y") || !plyIn.hasProperty("z")) {
    //     std::cerr << "PLY file does not contain vertex data!" << std::endl;
    //     return false;
    // }

    // Extract vertex positions
    std::vector<float> x = plyIn.getElement("vertex").getProperty<float>("x");
    std::vector<float> y = plyIn.getElement("vertex").getProperty<float>("y");
    std::vector<float> z = plyIn.getElement("vertex").getProperty<float>("z");

    // Resize and fill the vertex matrix V
    V.resize(x.size(), 3);
    for (size_t i = 0; i < x.size(); i++) {
        V.row(i) << x[i], y[i], z[i];
    }

    // Extract vertex normals (if available)
    if (plyIn.getElement("vertex").hasProperty("nx") &&
        plyIn.getElement("vertex").hasProperty("ny") &&
        plyIn.getElement("vertex").hasProperty("nz")) {
        std::vector<float> nx = plyIn.getElement("vertex").getProperty<float>("nx");
        std::vector<float> ny = plyIn.getElement("vertex").getProperty<float>("ny");
        std::vector<float> nz = plyIn.getElement("vertex").getProperty<float>("nz");

        // Resize and fill the normal matrix N
        N.resize(nx.size(), 3);
        for (size_t i = 0; i < nx.size(); i++) {
            N.row(i) << nx[i], ny[i], nz[i];
        }
    } else {
        std::cerr << "PLY file does not contain normal data!" << std::endl;
        N.resize(0, 3); // Clear N if no normals are present
    }

    return true;
}

void save_combined_spheres(
    const std::vector<Sphere>& circumspheres, // List of circumspheres (center, radius)
    const std::string& output_file)
{
    int resolution = 10;  // Sphere resolutiosn
    Eigen::MatrixXd V_base;
    Eigen::MatrixXi F_base;
    generate_sphere_mesh(1.0, resolution, V_base, F_base); // Generate unit sphere

    std::vector<Eigen::MatrixXd> vertices_list;
    std::vector<Eigen::MatrixXi> faces_list;
    int vertex_offset = 0;

    // Iterate through all circumspheres
    for (size_t i = 0; i < circumspheres.size(); ++i) {
        const cgp::vec3& p = circumspheres[i].center;
        double radius = circumspheres[i].radius;
		Point center(p.x, p.y, p.z);

        // Scale and translate the base sphere
        Eigen::MatrixXd V_scaled = V_base * radius; 
        Eigen::MatrixXd V_translated = V_scaled.rowwise() + Eigen::RowVector3d(center.x(), center.y(), center.z());

        // Adjust face indices
        Eigen::MatrixXi F_translated = F_base.array() + vertex_offset;

        // viewer.append_mesh();
        // viewer.data().set_mesh(V_translated, F_base);
        // viewer.data().set_colors(Eigen::RowVector3d(1,1,0));

        // Store the results
        vertices_list.push_back(V_translated);
        faces_list.push_back(F_translated);
        vertex_offset += V_translated.rows();
    }


    // Combine all vertices and faces into final matrices
    int total_vertices = 0;
    int total_faces = 0;

    for (size_t i = 0; i < vertices_list.size(); ++i) {
        total_vertices += vertices_list[i].rows();
        total_faces += faces_list[i].rows();
    }

    Eigen::MatrixXd V_final(total_vertices, 3);
    Eigen::MatrixXi F_final(total_faces, 3);

    int v_offset = 0, f_offset = 0;
    for (size_t i = 0; i < vertices_list.size(); ++i) {
        V_final.block(v_offset, 0, vertices_list[i].rows(), 3) = vertices_list[i];
        F_final.block(f_offset, 0, faces_list[i].rows(), 3) = faces_list[i];
        v_offset += vertices_list[i].rows();
        f_offset += faces_list[i].rows();
    }

    // Save the final mesh to an OBJ file
    igl::writeOBJ(output_file, V_final, F_final);
    std::cout << "Saved " << circumspheres.size() << " spheres to " << output_file << std::endl;
}


// Function to compute the circumsphere of a tetrahedron
Sphere compute_circumsphere(const Tetrahedron& tetra) {
    Eigen::Matrix3d A;
    Eigen::Vector3d b;

    Eigen::Vector3d p0(tetra.vertex(0).x(), tetra.vertex(0).y(), tetra.vertex(0).z());
    for (int i = 1; i < 4; ++i) {
        Eigen::Vector3d pi(tetra.vertex(i).x(), tetra.vertex(i).y(), tetra.vertex(i).z());
        A.row(i - 1) = 2 * (pi - p0);
        b(i - 1) = pi.squaredNorm() - p0.squaredNorm();
    }

    Eigen::Vector3d circumcenter = A.colPivHouseholderQr().solve(b);
    double radius = (circumcenter - p0).norm();
    return {cgp::vec3(circumcenter.x(), circumcenter.y(), circumcenter.z()), radius};
}

// Improved function to check if a point is inside the convex hull
bool is_inside_convex_hull(const Eigen::MatrixXd& V_hull,
                           const Eigen::MatrixXi& F_hull,
                           const Point& point) {
    for (int i = 0; i < F_hull.rows(); ++i) {
        Eigen::Vector3d v0 = V_hull.row(F_hull(i, 0));
        Eigen::Vector3d v1 = V_hull.row(F_hull(i, 1));
        Eigen::Vector3d v2 = V_hull.row(F_hull(i, 2));

        Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0).normalized();
        Eigen::Vector3d point_eigen(point.x(), point.y(), point.z());
        double signed_distance = normal.dot(point_eigen - v0);

        if (signed_distance > 1e-6) {
            return false;
        }
    }
    return true;
}

// Function to compute the angle between two vectors
double compute_angle(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    double dot_product = v1.dot(v2);
    return std::acos(dot_product / (v1.norm() * v2.norm()));
}

// / Function to process the mesh from the filepath and filter valid circumspheres
void scene_structure::filter_valid_circumspheres() {
    // Load the mesh
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
	Eigen::MatrixXd N;
	std::string ply_file=file_name;
	std::string normals_file;
	std::string output_file;
	// Automatically generate the normals filename by inserting "_normals" before ".ply"
    size_t dot_pos = ply_file.find_last_of(".");
    if (dot_pos != std::string::npos && ply_file.substr(dot_pos) == ".ply") {
        normals_file = ply_file.substr(0, dot_pos) + "_normals.ply";
        output_file =  ply_file.substr(0, dot_pos) + "_circumspheres.obj";

    } else {
        std::cerr << "Invalid PLY filename: " << ply_file << std::endl;
    }

    if (!igl::readPLY(ply_file, V, F)) {
        std::cerr << "Couldn't read PLY file: " << ply_file << std::endl;
    }

    if (!igl::readPLY(normals_file, N, F)) {
        std::cerr << "Couldn't read PLY file: " << normals_file << std::endl;
        
        readPLYWithNormals( ply_file,V, N) ;

    }

    std::cout << "Read " << V.rows() << " vertices and " << N.rows() << " Normals" << std::endl;
    
    assert(V.rows() == N.rows() && "Number of vertices and normals must be the same");

    igl::readPLY(file_name, V, F);
	std::vector<Point> points;
	std::cout << "Loaded mesh with " << V.rows() << " vertices and " << F.rows() << " faces" << std::endl;
    // Compute convex hull
    Eigen::MatrixXd V_hull;
    Eigen::MatrixXi F_hull;
    igl::copyleft::cgal::convex_hull(V, V_hull, F_hull);
	std::unordered_set<int> used_vertices;
	std::unordered_map<int, std::pair<Tetrahedron, double>> unused_vertices_info; // Maps unused vertex index to its original tetrahedron and dot product
	std::vector<Sphere> valid_spheres;

    // Prepare Delaunay triangulation
    Delaunay dt;
    for (int i = 0; i < V.rows(); ++i) {
        Point p(V(i, 0), V(i, 1), V(i, 2));
        points.push_back(p);
        dt.insert(p);
    }

    // Store tetrahedra and vertex indices
    std::vector<Tetrahedron> tetrahedra;
    std::vector<std::array<int, 4>> vertex_indices;
    std::vector<Sphere> circumspheres;
    for (auto it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it) {
        Tetrahedron tetra(it->vertex(0)->point(), it->vertex(1)->point(),
                        it->vertex(2)->point(), it->vertex(3)->point());

        // Get indices of the vertices of the current tetrahedron
        int v0 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(0)->point()));
        int v1 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(1)->point()));
        int v2 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(2)->point()));
        int v3 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(3)->point()));

        std::array<int, 4> vertex_indices = {v0, v1, v2, v3};

        auto [circumcenter, circumradius] = compute_circumsphere(tetra);

        // Check if circumcenter is inside the convex hull
        if (!is_inside_convex_hull(V_hull, F_hull, Point(circumcenter.x, circumcenter.y, circumcenter.z))) {
            continue;  // Skip this tetrahedron
        }

        bool all_valid = true;
        std::unordered_map<int, double> vertex_dot_products;

        for (size_t i = 0; i < 4; i++) {
            int vertex_idx = vertex_indices[i];
            // std::cout << "Vertex index: " << vertex_idx << std::endl;
            Eigen::RowVector3d vertex = Eigen::RowVector3d(tetra.vertex(i).x(),
                                                  tetra.vertex(i).y(),
                                                  tetra.vertex(i).z());

            Eigen::RowVector3d center=Eigen::RowVector3d(circumcenter.x, circumcenter.y, circumcenter.z);
            // Compute dot product for the vertex
            // double dot_product = N.row(vertex_idx).dot(V.row(vertex_idx)-center);

            double angle = compute_angle(vertex-center, N.row(vertex_idx));

            // If the angle is greater than 90 degrees, we reject it (we need the normal pointing towards the circumcenter)
            vertex_dot_products[vertex_idx] = angle;
            if (angle > 0.9) {
                all_valid = false;
            }
        }

        if (all_valid) {
            circumspheres.push_back({circumcenter, circumradius});
			field_function.spheres.push_back({circumcenter, circumradius});
			// std::cout << "Valid sphere added" << std::endl;

            // Mark vertices as used
            for (int vertex_idx : vertex_indices) {
                used_vertices.insert(vertex_idx);
            }
        } else {
            // Store unused vertices and their associated tetrahedron + dot product
            for (int vertex_idx : vertex_indices) {
                if (used_vertices.find(vertex_idx) == used_vertices.end()) { // If not used
                    unused_vertices_info[vertex_idx] = {tetra, vertex_dot_products[vertex_idx]};
                }
            }
        }
    }

    // Collect unused vertices
    std::vector<int> unused_vertices;
    for (int i = 0; i < V.rows(); ++i) {
        if (used_vertices.find(i) == used_vertices.end()) {
            unused_vertices.push_back(i);
        }
    }
	save_combined_spheres(circumspheres, output_file);


    // return valid_spheres;
}


void scene_structure::initialize()
{
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.set_rotation_axis_y();
	camera_control.look_at({ 3.0f, 2.0f, 2.0f }, {0,0,0}, {0,0,1});
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());


	// Initialization for the Implicit Surface
	// ***************************************** //
	// set the domain of the implicit surface boxx
	std::cout << "Set the domain of the implicit surface" <<gui.domain.samples<<gui.domain.length <<std::endl;
	
	implicit_surface.set_domain(gui.domain.samples, gui.domain.length);
	filter_valid_circumspheres();
	// field_function.spheres = { {field_function.a, 1.0}, {field_function.b, 1.0}, {field_function.c, 1.0} };
	std::cout << "Circumspheres filtered" <<field_function.spheres.size()<< std::endl;
	implicit_surface.update_field(field_function, gui.isovalue);
}



void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	
	// The standard frame
	if (gui.display.frame)
		draw(global_frame, environment);

	if (gui.display.surface)    // Display the implicit surface
		draw(implicit_surface.drawable_param.shape, environment);

	if (gui.display.wireframe) // Display the wireframe of the implicit surface
		draw_wireframe(implicit_surface.drawable_param.shape, environment, { 0,0,0 });

	if (gui.display.domain)    // Display the boundary of the domain
		draw(implicit_surface.drawable_param.domain_box, environment);

}

void scene_structure::display_gui()
{
	// Handle the gui values and the updates using the helper methods (*)
	implicit_surface.gui_update(gui, field_function);
}

void scene_structure::mouse_move_event()
{
	if (!inputs.keyboard.shift)
		camera_control.action_mouse_move(environment.camera_view);
}
void scene_structure::mouse_click_event()
{
	camera_control.action_mouse_click(environment.camera_view);
}
void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
}

