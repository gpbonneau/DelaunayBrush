#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <cmath>

#include <sstream>
#include <vector>
#include <string>
#include "circumspheredll.h"


// Define CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Tetrahedron_3 Tetrahedron;

// Function to compute the circumsphere of a tetrahedron
// std::pair<Point, double> compute_circumsphere(const Tetrahedron& tetra) {
//     const Point& p0 = tetra.vertex(0);
//     const Point& p1 = tetra.vertex(1);
//     const Point& p2 = tetra.vertex(2);
//     const Point& p3 = tetra.vertex(3);

//     // Compute the circumcenter using CGAL's circumcenter function
//     Point circumcenter = CGAL::circumcenter(p0, p1, p2, p3);
//     double radius = std::sqrt(CGAL::squared_distance(circumcenter, p0));

//     return {circumcenter, radius};
// }

// // Function to compute the angle between two vectors using dot product
// double compute_angle(const Vector& v1, const Vector& v2) {
//     double dot_product = v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
//     double magnitude_v1 = std::sqrt(v1.x() * v1.x() + v1.y() * v1.y() + v1.z() * v1.z());
//     double magnitude_v2 = std::sqrt(v2.x() * v2.x() + v2.y() * v2.y() + v2.z() * v2.z());
//     return std::acos(dot_product / (magnitude_v1 * magnitude_v2));  // In radians
// }

// Function to check if a point is inside any Delaunay tetrahedron
bool is_inside_delaunay(const Delaunay& dt, const Point& point) {
    auto cell = dt.locate(point);
    if (dt.is_infinite(cell)) {
        return false;  // Point is outside the convex hull
    }
    return true;  // Point is inside a finite tetrahedron
}

// // Function to compute valid circumspheres
// std::vector<std::pair<Point, double>> compute_valid_circumspheres(
//     const std::vector<Point>& points,
//     const std::vector<Vector>& normals) {

//     std::vector<std::pair<Point, int>> points_with_info;
//     for (int i = 0; i < points.size(); ++i) {
//         points_with_info.emplace_back(points[i], i);  // Assign index as info
//     }
//     Delaunay dt;
//     dt.insert(points_with_info.begin(), points_with_info.end());

//     std::vector<std::pair<Point, double>> valid_spheres;
//     std::unordered_set<int> used_vertices;

//     for (auto it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it) {
//         Tetrahedron tetra(it->vertex(0)->point(), it->vertex(1)->point(),
//                          it->vertex(2)->point(), it->vertex(3)->point());

//         auto [circumcenter, circumradius] = computeCircumspheres(tetra);

//         // Check if circumcenter is inside any Delaunay tetrahedron
//         if (!is_inside_delaunay(dt, circumcenter)) {
//             continue;  // Skip this tetrahedron
//         }
//         std::cout << "Circumcenter is inside a tetrahedron" << std::endl;

//         bool all_valid = true;
//         for (int i = 0; i < 4; ++i) {
//             int vertex_idx = it->vertex(i)->info();
//             // Get the vertex as a Point
//             Point vertex_point = it->vertex(i)->point();

//             // Compute the vector from the circumcenter to the vertex
//             Vector to_v = vertex_point - circumcenter;

//             double angle = compute_angle(to_v, normals[vertex_idx]);
//             if (angle > 0.9) {  // Angle greater than 90 degrees
//                 all_valid = false;
//                 break;
//             }
//         }
//         std::cout << "All vertices are valid: " << all_valid << std::endl;

//         if (all_valid) {
//             valid_spheres.push_back({circumcenter, circumradius});
//             for (int i = 0; i < 4; ++i) {
//                 used_vertices.insert(it->vertex(i)->info());
//             }
//         }
//     }

//     return valid_spheres;
// }

// Function to generate a sphere mesh
void generate_sphere_mesh(double radius, const Point& center, int resolution,
    std::vector<Point>& vertices, std::vector<std::array<int, 3>>& faces) {
vertices.clear();
faces.clear();

// Generate vertices
for (int i = 0; i <= resolution; ++i) {
double theta = M_PI * i / resolution; // Angle from top to bottom
for (int j = 0; j <= resolution; ++j) {
double phi = 2 * M_PI * j / resolution; // Angle around the sphere
double x = radius * sin(theta) * cos(phi);
double y = radius * sin(theta) * sin(phi);
double z = radius * cos(theta);

// Translate to the sphere's center
vertices.emplace_back(center.x() + x, center.y() + y, center.z() + z);
}
}

// Generate faces
for (int i = 0; i < resolution; ++i) {
for (int j = 0; j < resolution; ++j) {
int v0 = i * (resolution + 1) + j;
int v1 = v0 + 1;
int v2 = v0 + resolution + 1;
int v3 = v2 + 1;

faces.push_back({v0, v1, v2});
faces.push_back({v2, v1, v3});
}
}
}

// Function to save spheres to an OBJ file
void save_spheres_to_obj(const std::vector<std::pair<Point, double>>& spheres, const std::string& output_file) {
std::ofstream out(output_file);
if (!out) {
std::cerr << "Failed to open file: " << output_file << std::endl;
return;
}

int resolution = 20; // Sphere resolution (number of subdivisions)
int vertex_offset = 0; // Track the offset for face indices

for (const auto& [center, radius] : spheres) {
std::vector<Point> vertices;
std::vector<std::array<int, 3>> faces;

// Generate sphere mesh
generate_sphere_mesh(radius, center, resolution, vertices, faces);

// Write vertices to OBJ file
for (const auto& vertex : vertices) {
out << "v " << CGAL::to_double(vertex.x()) << " "
 << CGAL::to_double(vertex.y()) << " "
 << CGAL::to_double(vertex.z()) << "\n";
}

// Write faces to OBJ file (adjust indices by vertex_offset)
for (const auto& face : faces) {
out << "f " << face[0] + 1 + vertex_offset << " "
 << face[1] + 1 + vertex_offset << " "
 << face[2] + 1 + vertex_offset << "\n";
}

// Update vertex offset for the next sphere
vertex_offset += vertices.size();
}

out.close();
}

// // Main function to compute and save circumspheres
// void compute_and_save_circumspheres(
//     const std::vector<Point>& points,
//     const std::vector<Vector>& normals,
//     const std::string& output_file) {

//     auto start = std::chrono::high_resolution_clock::now();

//     // Compute valid circumspheres
//     auto valid_spheres = compute_valid_circumspheres(points, normals);

//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double, std::milli> duration = end - start;

//     std::cout << "Found " << valid_spheres.size() << " valid circumspheres\n";
//     std::cout << "Execution time: " << duration.count() << " ms\n";

//     // Save spheres to OBJ file
//     save_spheres_to_obj(valid_spheres, output_file);
//     std::cout << "Saved circumspheres to " << output_file << "\n";
// }



// Function to read vertices and normals from an OBJ file
bool read_obj(const std::string& filename, std::vector<Point>& points, std::vector<Vector>& normals) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "v") {
            // Read vertex
            double x, y, z;
            iss >> x >> y >> z;
            points.emplace_back(x, y, z);
        } else if (type == "vn") {
            // Read normal
            double nx, ny, nz;
            iss >> nx >> ny >> nz;
            normals.emplace_back(nx, ny, nz);
        }
    }

    file.close();

    if (points.size() != normals.size()) {
        std::cerr << "Mismatch between the number of vertices and normals in the OBJ file." << std::endl;
        return false;
    }

    return true;
}
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.obj>" << std::endl;
        return 1;
    }

    std::vector<Point> points;
    std::vector<Vector> normals;

    std::string obj_file = argv[1];
    if (!read_obj(obj_file, points, normals)) {
        std::cerr << "Failed to read OBJ file: " << obj_file << std::endl;
        return 1;
    }

    std::cout << "Read " << points.size() << " vertices and normals " << normals.size() << std::endl;

    // Convert to Point3D format for DLL
    std::vector<Point3D> input_points(points.size());
    std::vector<Point3D> input_normals(normals.size());

    for (size_t i = 0; i < points.size(); ++i) {
        input_points[i] = {CGAL::to_double(points[i].x()), CGAL::to_double(points[i].y()), CGAL::to_double(points[i].z())};
        input_normals[i] = {CGAL::to_double(normals[i].x()), CGAL::to_double(normals[i].y()), CGAL::to_double(normals[i].z())};
    }

    // Allocate output arrays
    int maxOutputSize = 10000;  // Adjust as needed
    std::vector<Point3D> outCenters(maxOutputSize);
    std::vector<double> outRadii(maxOutputSize);

    // Call the DLL function
    int num_spheres = computeCircumspheres(
        input_points.data(), input_normals.data(),
        static_cast<int>(input_points.size()),
        outCenters.data(), outRadii.data(), maxOutputSize);

    std::cout << "Computed " << num_spheres << " valid circumspheres\n";

    // Store results in a vector
    std::vector<std::pair<Point, double>> spheres;
    for (int i = 0; i < num_spheres; ++i) {
        spheres.emplace_back(Point(outCenters[i].x, outCenters[i].y, outCenters[i].z), outRadii[i]);
    }

    // Save results
    std::string output_file = obj_file.substr(0, obj_file.find_last_of(".")) + "_circumspheres.obj";
    save_spheres_to_obj(spheres, output_file);
    std::cout << "Saved circumspheres to " << output_file << "\n";

    return 0;
}
