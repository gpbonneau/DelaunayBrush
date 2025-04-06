


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <igl/harmonic.h>

#include <sstream>
#include <vector>
#include <string>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include "circumspheredll.h"

#include <igl/copyleft/marching_cubes.h>
#include <Eigen/Dense>
#include <vector>
#include <limits>

#include <igl/writeOBJ.h>   
#include <igl/readPLY.h>

#include <set>
#include <stdexcept>

#include <igl/readOBJ.h>
#include <igl/biharmonic_coordinates.h>
#include <igl/point_mesh_squared_distance.h>


// typedef CGAL::Nef_polyhedron_3<K> Nef_polyhedron;

typedef std::tuple<int, int, int> Voxel;
struct VoxelHash {
    std::size_t operator()(const Voxel& v) const {
        auto [x, y, z] = v;
        return std::hash<int>()(x) ^ (std::hash<int>()(y) << 1) ^ (std::hash<int>()(z) << 2);
    }
};





// Define CGAL types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Tetrahedron_3 Tetrahedron;
typedef CGAL::Surface_mesh<Point> Mesh;
struct PointHash {
    std::size_t operator()(const Point& p) const {
        std::size_t h1 = std::hash<double>()(CGAL::to_double(p.x()));
        std::size_t h2 = std::hash<double>()(CGAL::to_double(p.y()));
        std::size_t h3 = std::hash<double>()(CGAL::to_double(p.z()));
        return h1 ^ (h2 << 1) ^ (h3 << 2);  // Combine the hash values for x, y, and z
    }
};

// Convert Point3D to Eigen::Vector3d
Eigen::Vector3d toEigen(const Point3D &p) {
    return Eigen::Vector3d(p.x, p.y, p.z);
}

// Define the implicit function for the union of spheres
double union_of_spheres(const Eigen::Vector3d &x, const std::vector<Point3D> &centers, const std::vector<double> &radii) {
    double min_val = std::numeric_limits<double>::max();
    for (size_t i = 0; i < centers.size(); ++i) {
        double val = (x - toEigen(centers[i])).squaredNorm() - radii[i] * radii[i];
        min_val = std::min(min_val, val);
    }
    return min_val;
}



// Function to voxelize the spheres
void voxelize_spheres(const std::vector<std::pair<Point, double>>& spheres, double voxel_size,
    std::unordered_set<Voxel, VoxelHash>& voxel_grid) {
    for (const auto& [center, radius] : spheres) {
    int min_x = std::floor((center.x() - radius) / voxel_size);
    int max_x = std::ceil((center.x() + radius) / voxel_size);
    int min_y = std::floor((center.y() - radius) / voxel_size);
    int max_y = std::ceil((center.y() + radius) / voxel_size);
    int min_z = std::floor((center.z() - radius) / voxel_size);
    int max_z = std::ceil((center.z() + radius) / voxel_size);

    for (int x = min_x; x <= max_x; ++x) {
    for (int y = min_y; y <= max_y; ++y) {
    for (int z = min_z; z <= max_z; ++z) {
    Point voxel_center(x * voxel_size, y * voxel_size, z * voxel_size);
    if (CGAL::squared_distance(voxel_center, center) <= radius * radius) {
        voxel_grid.insert({x, y, z});
    }
    }
    }
    }
    }
    }

    // Function to extract boundary faces and save as OBJ
    void extract_and_save_boundary_faces(const std::unordered_set<Voxel, VoxelHash>& voxel_grid,
        double voxel_size, const std::string& filename) {
    std::vector<std::array<int, 3>> offsets = {
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}
    };

    std::vector<std::array<Point, 3>> faces;
    std::unordered_map<Point, int, PointHash> point_index;
    int index = 1;

    std::ofstream out(filename);
    if (!out) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
    }

    for (const auto& voxel : voxel_grid) {
    auto [x, y, z] = voxel;
    for (size_t i = 0; i < offsets.size(); ++i) {
    auto& offset = offsets[i];
    Voxel neighbor = {x + offset[0], y + offset[1], z + offset[2]};
    if (voxel_grid.find(neighbor) == voxel_grid.end()) {
    Point base(x * voxel_size, y * voxel_size, z * voxel_size);
    double half = voxel_size / 2.0;

    std::array<Point, 8> corners = {
    Point(base.x() - half, base.y() - half, base.z() - half),
    Point(base.x() + half, base.y() - half, base.z() - half),
    Point(base.x() + half, base.y() + half, base.z() - half),
    Point(base.x() - half, base.y() + half, base.z() - half),
    Point(base.x() - half, base.y() - half, base.z() + half),
    Point(base.x() + half, base.y() - half, base.z() + half),
    Point(base.x() + half, base.y() + half, base.z() + half),
    Point(base.x() - half, base.y() + half, base.z() + half)
    };

    std::array<std::array<int, 3>, 2> faces_indices;
    switch (i) {
    case 0: faces_indices = {{{0, 4, 7}, {7, 3, 0}}}; break;
    case 1: faces_indices = {{{1, 2, 6}, {6, 5, 1}}}; break;
    case 2: faces_indices = {{{0, 1, 5}, {5, 4, 0}}}; break;
    case 3: faces_indices = {{{2, 3, 7}, {7, 6, 2}}}; break;
    case 4: faces_indices = {{{0, 3, 2}, {2, 1, 0}}}; break;
    case 5: faces_indices = {{{4, 5, 6}, {6, 7, 4}}}; break;
    }

    for (const auto& indices : faces_indices) {
    faces.push_back({corners[indices[0]], corners[indices[1]], corners[indices[2]]});
    }
    }
    }
    }

    for (const auto& face : faces) {
    for (const auto& vertex : face) {
    if (point_index.find(vertex) == point_index.end()) {
    point_index[vertex] = index++;
    out << "v " << CGAL::to_double(vertex.x()) << " "
    << CGAL::to_double(vertex.y()) << " "
    << CGAL::to_double(vertex.z()) << "\n";
    }
    }
    }

    for (const auto& face : faces) {
    out << "f " << point_index[face[0]] << " "
    << point_index[face[1]] << " "
    << point_index[face[2]] << "\n";
    }

    out.close();
    std::cout << "Saved boundary voxelized mesh to " << filename << "\n";
}
// Function to generate a single sphere mesh

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
            normals.emplace_back(-1*nx, -1*ny, -1*nz);
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
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input.obj>" << std::endl;
        return 1;
    }

    std::vector<Point> points;
    std::vector<Vector> normals;

    std::string obj_file = argv[1];
    int res=std::stoi(argv[2]);


    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    // Save results
    std::string output_file = obj_file.substr(0, obj_file.find_last_of(".")) + "_circumspheres.obj";
    std::string filtered = obj_file.substr(0, obj_file.find_last_of(".")) + "_filtered.obj";
    std::string deformed = obj_file.substr(0, obj_file.find_last_of(".")) + "_final.obj";

    if(res>10)
    {
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

        
        // save_spheres_to_obj(spheres, output_file);

        
        Eigen::VectorXd S(res * res * res);
        Eigen::MatrixXd GV(res * res * res, 3);
        float mesh_size = 6.0;
        // Compute scalar field
        int index = 0;
        for (int i = 0; i < res; ++i) {
            for (int j = 0; j < res; ++j) {
                for (int k = 0; k < res; ++k) {
                    Eigen::Vector3d p(i * mesh_size / res - 1, j * mesh_size/ res - 1, k * mesh_size/ res - 1);
                    GV.row(index) = p;
                    S(index) = union_of_spheres(p, outCenters, outRadii);
                    index++;
                }
            }
        }

        // Extract mesh using marching cubes
        igl::copyleft::marching_cubes(S, GV, res, res, res, V, F);

        // Save the result
        igl::writeOBJ(output_file, V, F);

        std::cout << "Saved circumspheres to " << output_file << "\n";

}

    Eigen::MatrixXd filte_v;
    Eigen::MatrixXi filter_f;
    if(res<10)
    {   std::cout << output_file << std::endl;
        igl::readOBJ(output_file, V, F);
        igl::readOBJ(obj_file, filte_v, filter_f);
        deformed = obj_file.substr(0, obj_file.find_last_of(".")) + "_final2.obj";
    }
    else
    {
        igl::readPLY("circumsphere_vertices.ply", filte_v, filter_f);
    }

    Eigen::VectorXi I;
    Eigen::MatrixXd C, sqrD;
    std::cout << "Closest vertex: " << F.rows() << std::endl;
    igl::point_mesh_squared_distance(filte_v, V, F, sqrD, I, C);

    std::set<int> used_indices;
    std::vector<int> unique_indices;
    std::vector<Eigen::RowVector3d> unique_positions;

    for (int i = 0; i < filte_v.rows(); ++i) {
        int face_idx = I(i);
        Eigen::RowVector3d point = filte_v.row(i);

        // Get vertex indices of the closest triangle
        int v0 = F(face_idx, 0);
        int v1 = F(face_idx, 1);
        int v2 = F(face_idx, 2);

        // Find closest vertex among the 3
        std::vector<int> verts = {v0, v1, v2};
        int closest_vid = -1;
        double min_dist = std::numeric_limits<double>::max();

        for (int vid : verts) {
            double dist = (point - V.row(vid)).squaredNorm();
            if (dist < min_dist) {
                min_dist = dist;
                closest_vid = vid;
                
            }

        }

        // Skip if already used
        if (used_indices.count(closest_vid)) continue;

        used_indices.insert(closest_vid);
        unique_indices.push_back(closest_vid);
        unique_positions.push_back(point);
    }
    std::cout << "Unique indices: " << unique_indices.size() << std::endl;

    // Convert to Eigen matrices
    Eigen::VectorXi b(unique_indices.size());
    Eigen::MatrixXd bc(unique_positions.size(), 3);
    std::cout << "Unique indices: " << unique_indices.size() << std::endl;
    for (int i = 0; i < unique_indices.size(); ++i) {
        b(i) = unique_indices[i];
        bc.row(i) = V.row(unique_indices[i]);
    }

    // Compute weights
    Eigen::MatrixXd W;
    igl::harmonic(V, F, b, bc, 2, W);


    // Compute new vertex positions
    Eigen::MatrixXd V_new = W ;

    // Step 5: Save the deformed mesh
    if (!igl::writeOBJ(deformed, V_new, F)) {
        std::cerr << "Failed to save deformed mesh!" << std::endl;
        return -1;
    }

    std::cout << "Deformation complete! Saved as 'deformed_mesh.obj'" << std::endl;

    return 0;
}
