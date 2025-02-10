// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Delaunay_triangulation_3.h>
// #include <CGAL/Tetrahedron_3.h>
// #include <Eigen/Dense>
// #include <vector>
// #include <iostream>
// #include <igl/readPLY.h>
// #include <igl/writePLY.h>
// #include <igl/copyleft/cgal/convex_hull.h>
// #include <igl/opengl/glfw/Viewer.h>


// // Define CGAL types
// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL::Delaunay_triangulation_3<K> Delaunay;
// typedef K::Point_3 Point;
// typedef K::Tetrahedron_3 Tetrahedron;

// // Function to compute the circumsphere of a tetrahedron
// std::pair<Point, double> compute_circumsphere(const Tetrahedron& tetra) {
//     Eigen::Matrix3d A;
//     Eigen::Vector3d b;

//     Eigen::Vector3d p0(tetra.vertex(0).x(), tetra.vertex(0).y(), tetra.vertex(0).z());
//     for (int i = 1; i < 4; ++i) {
//         Eigen::Vector3d pi(tetra.vertex(i).x(), tetra.vertex(i).y(), tetra.vertex(i).z());
//         A.row(i - 1) = 2 * (pi - p0);
//         b(i - 1) = pi.squaredNorm() - p0.squaredNorm();
//     }

//     Eigen::Vector3d circumcenter = A.colPivHouseholderQr().solve(b);
//     double radius = (circumcenter - p0).norm();
//     return {Point(circumcenter.x(), circumcenter.y(), circumcenter.z()), radius};
// }

// // Function to check if a point is inside the convex hull
// bool is_inside_convex_hull(const Eigen::MatrixXd& V_hull,
//                            const Eigen::MatrixXi& F_hull,
//                            const Point& point) {
//     for (int i = 0; i < F_hull.rows(); ++i) {
//         // Get the vertices of the current face
//         Eigen::Vector3d v0 = V_hull.row(F_hull(i, 0));
//         Eigen::Vector3d v1 = V_hull.row(F_hull(i, 1));
//         Eigen::Vector3d v2 = V_hull.row(F_hull(i, 2));

//         // Compute the normal of the plane formed by the face
//         Eigen::Vector3d normal = (v1 - v0).cross(v2 - v0);
//         normal.normalize();

//         // Compute the signed distance of the point to the plane
//         Eigen::Vector3d point_eigen(point.x(), point.y(), point.z());
//         double signed_distance = normal.dot(point_eigen - v0);

//         // If the signed distance is positive, the point is outside the face
//         if (signed_distance > 1e-6) {
//             return false;  // The point is outside the convex hull
//         }
//     }

//     // If the point satisfies all face checks, it is inside the convex hull
//     return true;
// }

// // Function to compute the angle between two vectors using dot product
// double compute_angle(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
//     double dot_product = v1.dot(v2);
//     double magnitude_v1 = v1.norm();
//     double magnitude_v2 = v2.norm();
//     return std::acos(dot_product / (magnitude_v1 * magnitude_v2));  // In radians
// }


// bool is_valid_circumsphere(const Eigen::MatrixXd& V,
//                             const Eigen::MatrixXi& F,
//                             const Eigen::RowVector3d& n,  // Normals of the point cloud
//                             const Point& circumcenter,
//                             const Tetrahedron& tetra) {
//     Eigen::Vector3d circumcenter_eigen(circumcenter.x(), circumcenter.y(), circumcenter.z());
//     bool is_valid = true;

//     // Check the angle with respect to each vertex of the tetrahedron
//     for (int i = 0; i < 4; ++i) {
//         Eigen::Vector3d vertex = Eigen::Vector3d(tetra.vertex(i).x(),
//                                                   tetra.vertex(i).y(),
//                                                   tetra.vertex(i).z());

//         // Compute the vector from the circumcenter to the vertex
//         Eigen::Vector3d circum_to_vertex = vertex - circumcenter_eigen;

//         // Retrieve the normal of the vertex from the normals matrix using the index in V
//         Eigen::Vector3d normal = n;  // Get the normal of the vertex by its index

//         // Compute the angle between the circumcenter-to-vertex vector and the vertex normal
//         double angle = compute_angle(circum_to_vertex, normal);

//         // If the angle is greater than 90 degrees, we reject it (we need the normal pointing towards the circumcenter)
//         if (angle > 0.5) {
//             is_valid = false;
//             break;
//         }
//     }

//     return is_valid;
// }


// void generate_sphere_mesh(
//     double radius, int resolution,
//     Eigen::MatrixXd& V, Eigen::MatrixXi& F)
// {
//     std::vector<Eigen::RowVector3d> vertices;
//     std::vector<Eigen::RowVector3i> faces;

//     // Create vertices
//     for (int i = 0; i <= resolution; ++i)
//     {
//         double theta = M_PI * i / resolution; // Angle from top to bottom
//         for (int j = 0; j <= resolution; ++j)
//         {
//             double phi = 2 * M_PI * j / resolution; // Angle around the sphere
//             double x = radius * sin(theta) * cos(phi);
//             double y = radius * sin(theta) * sin(phi);
//             double z = radius * cos(theta);
//             vertices.emplace_back(x, y, z);
//         }
//     }

//     // Create faces
//     for (int i = 0; i < resolution; ++i)
//     {
//         for (int j = 0; j < resolution; ++j)
//         {
//             int v0 = i * (resolution + 1) + j;
//             int v1 = v0 + 1;
//             int v2 = v0 + resolution + 1;
//             int v3 = v2 + 1;

//             faces.emplace_back(v2, v1, v0);
//             faces.emplace_back(v2, v3, v1);
//         }
//     }

//     // Convert to Eigen matrices
//     V.resize(vertices.size(), 3);
//     F.resize(faces.size(), 3);

//     for (size_t i = 0; i < vertices.size(); ++i)
//         V.row(i) = vertices[i];

//     for (size_t i = 0; i < faces.size(); ++i)
//         F.row(i) = faces[i];
// }

// void save_combined_spheres(
//     const std::vector<std::pair<Point, double>>& circumspheres, // List of circumspheres (center, radius)
//     igl::opengl::glfw::Viewer& viewer,
//     const std::string& output_file)
// {
//     int resolution = 20;  // Sphere resolution
//     Eigen::MatrixXd V_base;
//     Eigen::MatrixXi F_base;
//     generate_sphere_mesh(1.0, resolution, V_base, F_base); // Generate unit sphere

//     std::vector<Eigen::MatrixXd> vertices_list;
//     std::vector<Eigen::MatrixXi> faces_list;
//     int vertex_offset = 0;

//     // Iterate through all circumspheres
//     for (size_t i = 0; i < circumspheres.size(); ++i) {
//         const Point& center = circumspheres[i].first;
//         double radius = circumspheres[i].second;

//         // Scale and translate the base sphere
//         Eigen::MatrixXd V_scaled = V_base * radius; 
//         Eigen::MatrixXd V_translated = V_scaled.rowwise() + Eigen::RowVector3d(center.x(), center.y(), center.z());

//         // Adjust face indices
//         Eigen::MatrixXi F_translated = F_base.array() + vertex_offset;

//         viewer.append_mesh();
//         viewer.data().set_mesh(V_translated, F_base);
//         viewer.data().set_colors(Eigen::RowVector3d(1,1,0));

//         // Store the results
//         vertices_list.push_back(V_translated);
//         faces_list.push_back(F_translated);
//         vertex_offset += V_translated.rows();
//     }


//     // Combine all vertices and faces into final matrices
//     int total_vertices = 0;
//     int total_faces = 0;

//     for (size_t i = 0; i < vertices_list.size(); ++i) {
//         total_vertices += vertices_list[i].rows();
//         total_faces += faces_list[i].rows();
//     }

//     Eigen::MatrixXd V_final(total_vertices, 3);
//     Eigen::MatrixXi F_final(total_faces, 3);

//     int v_offset = 0, f_offset = 0;
//     for (size_t i = 0; i < vertices_list.size(); ++i) {
//         V_final.block(v_offset, 0, vertices_list[i].rows(), 3) = vertices_list[i];
//         F_final.block(f_offset, 0, faces_list[i].rows(), 3) = faces_list[i];
//         v_offset += vertices_list[i].rows();
//         f_offset += faces_list[i].rows();
//     }

//     // Save the final mesh to an OBJ file
//     igl::writeOBJ(output_file, V_final, F_final);
// }


// int main(int argc, char** argv) {
//     if (argc < 2) {
//         std::cerr << "Usage: " << argv[0] << " <path_to_ply_file>" << std::endl;
//         return 1;
//     }

//    // Read the point cloud from the PLY file using libigl
//     Eigen::MatrixXd V; // Points (vertices)
//     Eigen::MatrixXi F; // Faces (not needed in this case, but required by libigl read function)
//     Eigen::MatrixXd N; // Normals (assumed to be in the PLY file)

//     std::string ply_file = argv[1];
//     std::string normals_file,output_file;

//     // Automatically generate the normals filename by inserting "_normals" before ".ply"
//     size_t dot_pos = ply_file.find_last_of(".");
//     if (dot_pos != std::string::npos && ply_file.substr(dot_pos) == ".ply") {
//         normals_file = ply_file.substr(0, dot_pos) + "_normals.ply";
//         output_file =  ply_file.substr(0, dot_pos) + "_circumspheres.obj";

//     } else {
//         std::cerr << "Invalid PLY filename: " << ply_file << std::endl;
//         return 1;
//     }

//     if (!igl::readPLY(ply_file, V, F)) {
//         std::cerr << "Couldn't read PLY file: " << ply_file << std::endl;
//         return 1;
//     }

//     if (!igl::readPLY(normals_file, N, F)) {
//         std::cerr << "Couldn't read PLY file: " << normals_file << std::endl;
//         return 1;
//     }

//     std::cout << "Read " << V.rows() << " vertices and " << N.rows() << " Normals" << std::endl;
    
//     assert(V.rows() == N.rows() && "Number of vertices and normals must be the same");
//     // Compute convex hull using libigl's CGAL wrapper
//     Eigen::MatrixXd V_hull;  // Convex hull vertices
//     Eigen::MatrixXi F_hull;  // Convex hull faces
//     igl::copyleft::cgal::convex_hull(V, V_hull, F_hull);

//     // Convert Eigen matrix to CGAL points for Delaunay triangulation
//     std::vector<Point> points;
//     for (int i = 0; i < V.rows(); ++i) {
//         points.push_back(Point(V(i, 0), V(i, 1), V(i, 2)));
//     }

//     // Compute Delaunay tetrahedralization
//     Delaunay dt;
//     dt.insert(points.begin(), points.end());
    
//     std::vector<std::pair<Point, double>> circumspheres;
//     std::unordered_map<int, std::pair<Tetrahedron, double>> unused_vertices_info; // Maps unused vertex index to its original tetrahedron and dot product
//     std::unordered_set<int> used_vertices; // Stores indices of vertices used in valid tetrahedra

//     for (auto it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it) {
//         Tetrahedron tetra(it->vertex(0)->point(), it->vertex(1)->point(),
//                         it->vertex(2)->point(), it->vertex(3)->point());

//         // Get indices of the vertices of the current tetrahedron
//         int v0 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(0)->point()));
//         int v1 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(1)->point()));
//         int v2 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(2)->point()));
//         int v3 = std::distance(points.begin(), std::find(points.begin(), points.end(), it->vertex(3)->point()));

//         std::array<int, 4> vertex_indices = {v0, v1, v2, v3};

//         auto [circumcenter, circumradius] = compute_circumsphere(tetra);

//         // Check if circumcenter is inside the convex hull
//         if (!is_inside_convex_hull(V_hull, F_hull, circumcenter)) {
//             continue;  // Skip this tetrahedron
//         }

//         bool all_valid = true;
//         std::unordered_map<int, double> vertex_dot_products;

//         for (size_t i = 0; i < 4; i++) {
//             int vertex_idx = vertex_indices[i];
//             // std::cout << "Vertex index: " << vertex_idx << std::endl;
//             Eigen::RowVector3d vertex = Eigen::RowVector3d(tetra.vertex(i).x(),
//                                                   tetra.vertex(i).y(),
//                                                   tetra.vertex(i).z());

//             Eigen::RowVector3d center=Eigen::RowVector3d(circumcenter.x(), circumcenter.y(), circumcenter.z());
//             // Compute dot product for the vertex
//             // double dot_product = N.row(vertex_idx).dot(V.row(vertex_idx)-center);

//             double angle = compute_angle(vertex-center, N.row(vertex_idx));

//             // If the angle is greater than 90 degrees, we reject it (we need the normal pointing towards the circumcenter)
//             vertex_dot_products[vertex_idx] = angle;
//             if (angle > 0.9) {
//                 all_valid = false;
//             }
//         }

//         if (all_valid) {
//             circumspheres.push_back({circumcenter, circumradius});

//             // Mark vertices as used
//             for (int vertex_idx : vertex_indices) {
//                 used_vertices.insert(vertex_idx);
//             }
//         } else {
//             // Store unused vertices and their associated tetrahedron + dot product
//             for (int vertex_idx : vertex_indices) {
//                 if (used_vertices.find(vertex_idx) == used_vertices.end()) { // If not used
//                     unused_vertices_info[vertex_idx] = {tetra, vertex_dot_products[vertex_idx]};
//                 }
//             }
//         }
//     }

//     // Collect unused vertices
//     std::vector<int> unused_vertices;
//     for (int i = 0; i < V.rows(); ++i) {
//         if (used_vertices.find(i) == used_vertices.end()) {
//             unused_vertices.push_back(i);
//         }
//     }

    
//     std::cout << "Total circumspheres stored: " << circumspheres.size() << std::endl;
    
//     // Prepare to write results to a PLY file
//     Eigen::MatrixXd circumsphere_points(circumspheres.size(), 3);
//     Eigen::MatrixXd radii(circumspheres.size(), 1);
    
//     igl::opengl::glfw::Viewer viewer;
// //     for (size_t i = 0; i < circumspheres.size(); ++i) {
// //         const Point& center = circumspheres[i].first;
// //         circumsphere_points(i, 0) = center.x();
// //         circumsphere_points(i, 1) = center.y();
// //         circumsphere_points(i, 2) = center.z();
// //         radii(i, 0) = circumspheres[i].second;

// //         Eigen::RowVector3d sphere_center = circumsphere_points.row(i);
// //         double radius = circumspheres[i].second;
// //         int resolution = 5;
// //         // Generate a base sphere mesh
// //         Eigen::MatrixXd V_base;
// //         Eigen::MatrixXi F_base;
// //         generate_sphere_mesh(radius, resolution, V_base, F_base);
// //         Eigen::MatrixXd V_translated = V_base.rowwise() + sphere_center;

// //         // Add the sphere to the viewer
// //         viewer.append_mesh();
// //         viewer.data().set_mesh(V_translated, F_base);
// //         viewer.data().set_colors(Eigen::RowVector3d(1,1,0));
// // }
//     std::cout << "Writing circumspheres to 'circumspheres.obj'" << std::endl;
//     save_combined_spheres(circumspheres,viewer, output_file);



//     // Write the circumsphere points and radii to a new PLY file
//     // igl::writePLY("circumspheres.ply", circumsphere_points, F);

//     std::cout << "Circumsphere points written to 'circumspheres.ply'" << std::endl;

//     viewer.data().add_points(V, Eigen::RowVector3d(1, 0, 0));
//     viewer.launch();
    
//     return 0;
// }



#include "cgp/cgp.hpp" // Give access to the complete CGP library
#include "environment.hpp" // The general scene environment + project variable
#include <iostream> 

#include <chrono>
#include <thread>

// Custom scene of this code
#include "scene.hpp"




// *************************** //
// Custom Scene defined in "scene.hpp"
// *************************** //

scene_structure scene;


// The rest of this code is a generic initialization and animation loop that can be applied to different scenes
// *************************** //
// Start of the program
// *************************** //

window_structure standard_window_initialization(int width = 0, int height = 0);
void initialize_default_shaders();
void animation_loop();
void display_gui_default();

timer_fps fps_record;

int main(int, char* argv[])
{
	
    std::string ply_file = argv[1];
    scene.file_name=ply_file;

	

	// ************************ //
	//     INITIALISATION
	// ************************ //

	// Standard Initialization of an OpenGL ready window
	scene.window = standard_window_initialization();

	// Initialize System Info
	project::path = cgp::project_path_find(argv[0], "shaders/");

	// Initialize default shaders
	initialize_default_shaders();


	// Custom scene initialization
	std::cout << "Initialize data of the scene ..." << std::endl;
	scene.initialize();
	std::cout << "Initialization finished\n" << std::endl;
	std::cout<<"\n\n"+scene.camera_control.doc_usage()<<std::endl;


	// ************************ //
	//     Animation Loop
	// ************************ //
	std::cout << "Start animation loop ..." << std::endl;
	fps_record.start();


	// Call the main display loop in the function animation_loop
	//  The following part is simply a loop that call the function "animation_loop"
	//  (This call is different when we compile in standard mode with GLFW, than when we compile with emscripten to output the result in a webpage.)
#ifndef __EMSCRIPTEN__
    double lasttime = glfwGetTime();
	// Default mode to run the animation/display loop with GLFW in C++
	while (!glfwWindowShouldClose(scene.window.glfw_window)) {
		// The real animation loop
		animation_loop();

		// FPS limitation
		if(project::fps_limiting){
			while (glfwGetTime() < lasttime + 1.0 / project::fps_max) {	}
        	lasttime = glfwGetTime();
		}
	}
#else
	// Specific loop if compiled for EMScripten
	emscripten_set_main_loop(animation_loop, 0, 1);
#endif

	std::cout << "\nAnimation loop stopped" << std::endl;

	// Cleanup
	cgp::imgui_cleanup();
	glfwDestroyWindow(scene.window.glfw_window);
	glfwTerminate();

	return 0;
}

void animation_loop()
{

	emscripten_update_window_size(scene.window.width, scene.window.height); // update window size in case of use of emscripten (not used by default)

	scene.camera_projection.aspect_ratio = scene.window.aspect_ratio();
	scene.environment.camera_projection = scene.camera_projection.matrix();
	glViewport(0, 0, scene.window.width, scene.window.height);

	vec3 const& background_color = scene.environment.background_color;
	glClearColor(background_color.x, background_color.y, background_color.z, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);

	float const time_interval = fps_record.update();
	if (fps_record.event) {
		std::string const title = "CGP Display - " + str(fps_record.fps) + " fps";
		glfwSetWindowTitle(scene.window.glfw_window, title.c_str());
	}

	imgui_create_frame();
	ImGui::GetIO().FontGlobalScale = project::gui_scale;
	ImGui::Begin("GUI", NULL, ImGuiWindowFlags_AlwaysAutoResize);
	scene.inputs.mouse.on_gui = ImGui::GetIO().WantCaptureMouse;
	scene.inputs.time_interval = time_interval;


	// Display the ImGUI interface (button, sliders, etc)
	display_gui_default();
	scene.display_gui();

	// Handle camera behavior in standard frame
	scene.idle_frame();

	// Call the display of the scene
	scene.display_frame();


	// End of ImGui display and handle GLFW events
	ImGui::End();
	imgui_render_frame(scene.window.glfw_window);
	glfwSwapBuffers(scene.window.glfw_window);
	glfwPollEvents();
}


void initialize_default_shaders()
{
	// Generate the default directory from which the shaders are found
	//  By default, it should be "shaders/"
	std::string default_path_shaders = project::path +"shaders/";

	// Set standard mesh shader for mesh_drawable
	mesh_drawable::default_shader.load(default_path_shaders +"mesh/mesh.vert.glsl", default_path_shaders +"mesh/mesh.frag.glsl");
	triangles_drawable::default_shader.load(default_path_shaders +"mesh/mesh.vert.glsl", default_path_shaders +"mesh/mesh.frag.glsl");

	// Set default white texture
	image_structure const white_image = image_structure{ 1,1,image_color_type::rgba,{255,255,255,255} };
	mesh_drawable::default_texture.initialize_texture_2d_on_gpu(white_image);
	triangles_drawable::default_texture.initialize_texture_2d_on_gpu(white_image);

	// Set standard uniform color for curve/segment_drawable
	curve_drawable::default_shader.load(default_path_shaders +"single_color/single_color.vert.glsl", default_path_shaders+"single_color/single_color.frag.glsl");
}





//Callback functions
void window_size_callback(GLFWwindow* window, int width, int height);
void mouse_move_callback(GLFWwindow* window, double xpos, double ypos);
void mouse_scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void mouse_click_callback(GLFWwindow* window, int button, int action, int mods);
void keyboard_callback(GLFWwindow* window, int key, int, int action, int mods);

// Standard initialization procedure
window_structure standard_window_initialization(int width_target, int height_target)
{
	// Create the window using GLFW
	// ***************************************************** //
	window_structure window;
	window.initialize(width_target, height_target, "CGP Display", CGP_OPENGL_VERSION_MAJOR, CGP_OPENGL_VERSION_MINOR);

	// Display information
	// ***************************************************** //

	// Display window size
	std::cout << "\nWindow (" << window.width << "px x " << window.height << "px) created" << std::endl;
	std::cout << "Monitor: " << glfwGetMonitorName(window.monitor) << " - Resolution (" << window.screen_resolution_width << "x" << window.screen_resolution_height << ")\n" << std::endl;

	// Display debug information on command line
	std::cout << "OpenGL Information:" << std::endl;
	std::cout << cgp::opengl_info_display() << std::endl;

	// Initialize ImGUI
	// ***************************************************** //
	cgp::imgui_init(window.glfw_window);

	// Set the callback functions for the inputs
	glfwSetMouseButtonCallback(window.glfw_window, mouse_click_callback); // Event called when a button of the mouse is clicked/released
	glfwSetCursorPosCallback(window.glfw_window, mouse_move_callback);    // Event called when the mouse is moved
	glfwSetWindowSizeCallback(window.glfw_window, window_size_callback);  // Event called when the window is rescaled        
	glfwSetKeyCallback(window.glfw_window, keyboard_callback);            // Event called when a keyboard touch is pressed/released
	glfwSetScrollCallback(window.glfw_window, mouse_scroll_callback);     // Event called when scrolling the mouse

	return window;
}




// This function is called everytime the window is resized
void window_size_callback(GLFWwindow*, int width, int height)
{
	scene.window.width = width;
	scene.window.height = height;
}

// This function is called everytime the mouse is moved
void mouse_move_callback(GLFWwindow* /*window*/, double xpos, double ypos)
{
	vec2 const pos_relative = scene.window.convert_pixel_to_relative_coordinates({ xpos, ypos });
	scene.inputs.mouse.position.update(pos_relative);
	scene.mouse_move_event();
}

// This function is called everytime a mouse button is clicked/released
void mouse_click_callback(GLFWwindow* window, int button, int action, int mods)
{
	ImGui_ImplGlfw_MouseButtonCallback(window, button, action, mods);
	
	scene.inputs.mouse.click.update_from_glfw_click(button, action);
	scene.mouse_click_event();
}

// This function is called everytime the mouse is scrolled
void mouse_scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	ImGui_ImplGlfw_ScrollCallback(window, xoffset, yoffset);

	scene.inputs.mouse.scroll = yoffset;
	scene.mouse_scroll_event();
}

// This function is called everytime a keyboard touch is pressed/released
void keyboard_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	ImGui_ImplGlfw_KeyCallback(window, key, scancode, action, mods);
	bool imgui_capture_keyboard = ImGui::GetIO().WantCaptureKeyboard;
	
	if(!imgui_capture_keyboard){
		scene.inputs.keyboard.update_from_glfw_key(key, action);
		scene.keyboard_event();

		// Press 'F' for full screen mode
		if (key == GLFW_KEY_F && action == GLFW_PRESS && scene.inputs.keyboard.shift) {
			scene.window.is_full_screen = !scene.window.is_full_screen;
			if (scene.window.is_full_screen)
				scene.window.set_full_screen();
			else
				scene.window.set_windowed_screen();
		}
		// Press 'V' for camera frame/view matrix debug
		if (key == GLFW_KEY_V && action == GLFW_PRESS && scene.inputs.keyboard.shift) {
			auto const camera_model = scene.camera_control.camera_model;
			std::cout << "\nDebug camera (position = " << str(camera_model.position()) << "):\n" << std::endl;
			std::cout << "  Frame matrix:" << std::endl;
			std::cout << str_pretty(camera_model.matrix_frame()) << std::endl;
			std::cout << "  View matrix:" << std::endl;
			std::cout << str_pretty(camera_model.matrix_view()) << std::endl;

		}
	}

}

void display_gui_default()
{
	std::string fps_txt = str(fps_record.fps)+" fps";
	ImGui::Text( fps_txt.c_str(), "%s" );
	if(ImGui::CollapsingHeader("Window")) {
#ifndef __EMSCRIPTEN__
		bool changed_screen_mode = ImGui::Checkbox("Full Screen", &scene.window.is_full_screen);
		if(changed_screen_mode){	
			if (scene.window.is_full_screen)
				scene.window.set_full_screen();
			else
				scene.window.set_windowed_screen();
		}
#endif
		ImGui::SliderFloat("Gui Scale", &project::gui_scale, 0.5f, 2.5f);

#ifndef __EMSCRIPTEN__
		// Arbitrary limits the refresh rate to a maximal frame per seconds.
		//  This limits the risk of having different behaviors when you use different machine. 
		ImGui::Checkbox("FPS limiting",&project::fps_limiting);
		if(project::fps_limiting){
			ImGui::SliderFloat("FPS limit",&project::fps_max, 10, 250);
		}
#endif
		// vsync is the default synchronization of frame refresh with the screen frequency
		//   vsync may or may not be enforced by your GPU driver and OS (on top of the GLFW request).
		//   de-activating vsync may generate arbitrary large FPS depending on your GPU and scene.
		if(ImGui::Checkbox("vsync (screen sync)",&project::vsync)){
			project::vsync==true? glfwSwapInterval(1) : glfwSwapInterval(0); 
		}


		ImGui::Spacing();ImGui::Separator();ImGui::Spacing();
	}
}


