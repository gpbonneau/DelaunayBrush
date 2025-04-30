#include "CircumsphereDLL.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <cmath>

// CGAL type definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;
typedef K::Tetrahedron_3 Tetrahedron;

// Compute circumsphere
std::pair<Point, double> compute_circumsphere(const Tetrahedron& tetra) {
    Point circumcenter = CGAL::circumcenter(
        tetra.vertex(0), tetra.vertex(1), tetra.vertex(2), tetra.vertex(3));
    double radius = std::sqrt(CGAL::squared_distance(circumcenter, tetra.vertex(0)));
    return {circumcenter, radius};
}

// Compute angle between two vectors
double compute_angle(const Vector& v1, const Vector& v2) {
    double dot_product = v1 * v2;
    double magnitude_v1 = std::sqrt(v1.squared_length());
    double magnitude_v2 = std::sqrt(v2.squared_length());
    return std::acos(dot_product / (magnitude_v1 * magnitude_v2));
}

// Export unique vertices to a point cloud file
void export_vertices_to_ply(const std::unordered_set<Point>& unique_vertices, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return;
    }
    
    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "element vertex " << unique_vertices.size() << "\n";
    file << "property float x\n";
    file << "property float y\n";
    file << "property float z\n";
    file << "end_header\n";
    
    for (const auto& point : unique_vertices) {
        file << point.x() << " " << point.y() << " " << point.z() << "\n";
    }
    
    file.close();
}

// Main DLL function to compute circumspheres
extern "C" DLL_EXPORT int computeCircumspheres(
    const Point3D* points, const Point3D* normals, int numPoints,
    Point3D* outCenters, double* outRadii, int maxOutputSize,std::string filename) {

    if (numPoints <= 0 || !points || !normals || !outCenters || !outRadii) {
        return 0;
    }

    std::vector<Point> cgal_points;
    std::vector<Vector> cgal_normals;
    std::unordered_set<Point> unique_vertices;

    for (int i = 0; i < numPoints; ++i) {
        Point p(points[i].x, points[i].y, points[i].z);
        cgal_points.push_back(p);
        cgal_normals.emplace_back(normals[i].x, normals[i].y, normals[i].z);
    }

    // Perform Delaunay triangulation
    Delaunay dt;
    std::vector<std::pair<Point, int>> points_with_info;
    for (int i = 0; i < numPoints; ++i) {
        points_with_info.emplace_back(cgal_points[i], i);
    }
    dt.insert(points_with_info.begin(), points_with_info.end());

    int count = 0;
    for (auto it = dt.finite_cells_begin(); it != dt.finite_cells_end(); ++it) {
        if (count >= maxOutputSize) break;

        Tetrahedron tetra(it->vertex(0)->point(), it->vertex(1)->point(),
                          it->vertex(2)->point(), it->vertex(3)->point());

        std::pair<Point, double> cir = compute_circumsphere(tetra);
        Point circumcenter = cir.first;
        double circumradius = cir.second;
        
        auto cell = dt.locate(circumcenter);
        if (dt.is_infinite(cell)) {
            continue;
        }

        bool all_valid = true;
        for (int i = 0; i < 4; ++i) {
            int vertex_idx = it->vertex(i)->info();
            Vector to_v = it->vertex(i)->point() - circumcenter;
            double angle = compute_angle(to_v, cgal_normals[vertex_idx]);
            if (angle > 0.9) {
                all_valid = false;
                break;
            }
        }

        if (all_valid) {
            outCenters[count] = {circumcenter.x(), circumcenter.y(), circumcenter.z()};
            outRadii[count] = circumradius;
            count++;
            for (int i = 0; i < 4; ++i) {
                unique_vertices.insert(it->vertex(i)->point());
            }
        }

        
    }
    
    export_vertices_to_ply(unique_vertices, filename);
    
    return count;
}
