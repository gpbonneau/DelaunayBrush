// #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #include <CGAL/Delaunay_triangulation_3.h>
// #include <CGAL/Triangulation_vertex_base_with_info_3.h>
// #include <CGAL/Tetrahedron_3.h>
// #include <Eigen/Dense>
// #include <vector>
// #include <iostream>
// #include <igl/readPLY.h>
// #include <igl/writePLY.h>
// #include <igl/copyleft/cgal/convex_hull.h>
// #include <igl/opengl/glfw/Viewer.h>
// #include <chrono>
// // #include <omp.h> // Add OpenMP header

// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
// typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;
// typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
// typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
// typedef K::Point_3 Point;
// typedef K::Plane_3 Plane;
// typedef K::Vector_3 Vector;
// typedef K::Tetrahedron_3 Tetrahedron;

// void generate_sphere_mesh(double radius, int resolution, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
//     V.resize((resolution+1)*(resolution+1), 3);
//     F.resize(2*resolution*resolution, 3);
//     for(int i=0; i<=resolution; ++i) {
//         double theta = M_PI * i / resolution;
//         for(int j=0; j<=resolution; ++j) {
//             double phi = 2 * M_PI * j / resolution;
//             int idx = i*(resolution+1) + j;
//             V.row(idx) << radius*sin(theta)*cos(phi), radius*sin(theta)*sin(phi), radius*cos(theta);
//         }
//     }
//     int f = 0;
//     for(int i=0; i<resolution; ++i) {
//         for(int j=0; j<resolution; ++j) {
//             int v0 = i*(resolution+1) + j;
//             int v1 = v0 + 1;
//             int v2 = v0 + (resolution+1);
//             int v3 = v2 + 1;
//             F.row(f++) << v0, v1, v2;
//             F.row(f++) << v2, v1, v3;
//         }
//     }
// }

// void save_combined_spheres(const std::vector<std::pair<Point, double>>& spheres, const std::string& output_file) {
//     Eigen::MatrixXd V_base;
//     Eigen::MatrixXi F_base;
//     generate_sphere_mesh(1.0, 10, V_base, F_base); // Lower resolution

//     Eigen::MatrixXd V_final;
//     Eigen::MatrixXi F_final;
//     int offset = 0;
//     for(const auto& [center, radius] : spheres) {
//         Eigen::MatrixXd V_scaled = V_base * radius;
//         V_scaled.rowwise() += Eigen::RowVector3d(CGAL::to_double(center.x()), 
//                                                 CGAL::to_double(center.y()), 
//                                                 CGAL::to_double(center.z()));
//         Eigen::MatrixXi F_offset = F_base.array() + offset;
//         V_final.conservativeResize(V_final.rows() + V_scaled.rows(), 3);
//         V_final.bottomRows(V_scaled.rows()) = V_scaled;
//         F_final.conservativeResize(F_final.rows() + F_offset.rows(), 3);
//         F_final.bottomRows(F_offset.rows()) = F_offset;
//         offset += V_scaled.rows();
//     }
//     igl::writeOBJ(output_file, V_final, F_final);
// }

// int main(int argc, char** argv) {
//     if(argc < 2) {
//         std::cerr << "Usage: " << argv[0] << " <input.ply>" << std::endl;
//         return 1;
//     }
//     auto start = std::chrono::high_resolution_clock::now();

//     Eigen::MatrixXd V, N;
//     Eigen::MatrixXi F;
//     std::string ply_file = argv[1];
//     if(!igl::readPLY(ply_file, V, F)) {
//         std::cerr << "Failed to read: " << ply_file << std::endl;
//         return 1;
//     }

//     size_t dot_pos = ply_file.find_last_of(".");
//     std::string normals_file = ply_file.substr(0, dot_pos) + "_normals.ply";
//     if(!igl::readPLY(normals_file, N, F)) {
//         std::cerr << "Failed to read: " << normals_file << std::endl;
//         return 1;
//     }

//     Eigen::MatrixXd V_hull;
//     Eigen::MatrixXi F_hull;
//     igl::copyleft::cgal::convex_hull(V, V_hull, F_hull);

//     std::vector<std::pair<Point, int>> points;
//     for(int i=0; i<V.rows(); ++i)
//         points.emplace_back(Point(V(i,0), V(i,1), V(i,2)), i);

//     Delaunay dt;
//     dt.insert(points.begin(), points.end());

//     std::vector<Plane> hull_planes;
//     for(int i=0; i<F_hull.rows(); ++i) {
//         Eigen::RowVector3d v0 = V_hull.row(F_hull(i,0));
//         Eigen::RowVector3d v1 = V_hull.row(F_hull(i,1));
//         Eigen::RowVector3d v2 = V_hull.row(F_hull(i,2));
//         hull_planes.emplace_back(Point(v0[0],v0[1],v0[2]), 
//                                 Point(v1[0],v1[1],v1[2]), 
//                                 Point(v2[0],v2[1],v2[2]));
//     }

//     std::vector<Vector> normals;
//     for(int i=0; i<N.rows(); ++i)
//         normals.emplace_back(N(i,0), N(i,1), N(i,2));

//     std::vector<std::pair<Point, double>> valid_spheres;
//     std::vector<std::vector<std::pair<Point, double>>> thread_spheres(omp_get_max_threads());

//     #pragma omp parallel
//     {
//         int thread_id = omp_get_thread_num();
//         #pragma omp for
//         for(auto cit = dt.finite_cells_begin(); cit != dt.finite_cells_end(); ++cit) {
//             Tetrahedron tetra(cit->vertex(0)->point(), cit->vertex(1)->point(),
//                              cit->vertex(2)->point(), cit->vertex(3)->point());
//             Point cc = CGAL::circumcenter(tetra);
            
//             bool inside = true;
//             for(const auto& plane : hull_planes) {
//                 if(plane.oriented_side(cc) != CGAL::ON_NEGATIVE_SIDE) {
//                     inside = false;
//                     break;
//                 }
//             }
//             if(!inside) continue;

//             bool valid = true;
//             for(int i=0; i<4; ++i) {
//                 int idx = cit->vertex(i)->info();
//                 Vector to_v = tetra.vertex(i) - cc;
//                 if(to_v * normals[idx] <= 0) {
//                     valid = false;
//                     break;
//                 }
//             }
//             if(valid) {
//                 double radius = CGAL::sqrt(CGAL::squared_distance(cc, tetra.vertex(0)));
//                 thread_spheres[thread_id].emplace_back(cc, radius);
//             }
//         }
//     }

//     // Merge results from all threads
//     for(const auto& vec : thread_spheres) {
//         valid_spheres.insert(valid_spheres.end(), vec.begin(), vec.end());
//     }

//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double, std::milli> duration = end - start;
//     std::cout << "Spheres found: " << valid_spheres.size() << "\nExecution time: " << duration.count() << " ms" << std::endl;

//     std::string output_file = ply_file.substr(0, dot_pos) + "_circumspheres.obj";
//     save_combined_spheres(valid_spheres, output_file);
//     std::cout << "Saved " << valid_spheres.size() << " spheres to " << output_file << std::endl;

//     return 0;
// }