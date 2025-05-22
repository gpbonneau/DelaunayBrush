#include <igl/readOBJ.h>   
#include <Eigen/Dense>
#include <Eigen/Sparse>


#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/polyscope.h"


#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"


// Track state statically outside of your callback (or bind like a lambda)
int myVal = 10;

void callSomeOtherFunction() {
    printf("coucou function!");
}

void myCallback() { // gets executed per-frame
    printf("coucou callback!");

  // Update content in the scene
//   polyscope::registerPointCloud("this frame point", my_points);

  // Build a UI element to edit a parameter, which will 
  // appear in the onscreen panel
  ImGui::InputInt("my val", &myVal); 

  if (ImGui::Button("run subroutine")) {
    callSomeOtherFunction();
  }
}

int main(int argc, char** argv) {
    Eigen::MatrixXd strokeV;
    Eigen::MatrixXi strokeF;

    std::string obj_file = argv[1];
    int N = atoi(argv[2]);

    // Save results
    std::string output_file = obj_file.substr(0, obj_file.find_last_of(".")) + "_marching.obj";
    std::string stroke_file = obj_file.substr(0, obj_file.find_last_of("_")) + "_strokes.obj";
    std::string dual_contour = obj_file.substr(0, obj_file.find_last_of(".")) + "_dual_contouring.obj";
    std::string sphere_file = obj_file.substr(0, obj_file.find_last_of(".")) + "_circumspheres.txt";
    std::string sphere_file_obj = obj_file.substr(0, obj_file.find_last_of(".")) + "_circumspheres.obj";
    std::string filtered = obj_file.substr(0, obj_file.find_last_of(".")) + "_filtered.ply";
    std::string final_mc = obj_file.substr(0, obj_file.find_last_of(".")) + "_final_mc.obj";
    std::string final_dc = obj_file.substr(0, obj_file.find_last_of(".")) + "_final_dc.obj";
    std::string timing = obj_file.substr(0, obj_file.find_last_of(".")) + "timing.txt";
    std::string uniorm_mesh = obj_file.substr(0, obj_file.find_last_of(".")) + "_uniform.obj";

    igl::readOBJ(stroke_file, strokeV, strokeF);

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;


    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, stroke_file.c_str());

    std::cout << "Loaded OBJ file: " << stroke_file << "\n";

    for (size_t i = 0; i < shapes.size(); ++i) {
        const tinyobj::shape_t& shape = shapes[i];
        std::cout << "Shape [" << i << "] Name: " << shape.name << "\n";
        std::cout << "  Number of faces: " << shape.mesh.num_face_vertices.size() << "\n";
    }

    // Load all vertex positions
    std::vector<std::array<double, 3>> vertices;
    for (size_t i = 0; i < attrib.vertices.size(); i += 3) {
        vertices.push_back({
            attrib.vertices[i + 0],
            attrib.vertices[i + 1],
            attrib.vertices[i + 2]
        });
    }


    // Collect faces from first N shapes
    std::vector<std::array<size_t, 3>> faces;
    for (int s = 0; s < std::min(N, (int)shapes.size()); ++s) {
        const auto& shape = shapes[s];
        size_t index_offset = 0;

        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); ++f) {
            int fv = shape.mesh.num_face_vertices[f];
            if (fv != 3) {
                index_offset += fv;
                continue; // only handle triangles
            }

            std::array<size_t, 3> face;
            for (size_t v = 0; v < 3; ++v) {
                face[v] = shape.mesh.indices[index_offset + v].vertex_index;
            }
            faces.push_back(face);
            index_offset += 3;
        }
    }

    // Initialize polyscope
    polyscope::init();
    // Register the mesh with Polyscope
    polyscope::view::setNavigateStyle(polyscope::NavigateStyle::Free);
    // polyscope::registerSurfaceMesh("Sculpted surface", mcV_new, mcF);
    // polyscope::registerSurfaceMesh("Strokes", strokeV, strokeF);
    polyscope::registerSurfaceMesh("Strokes", vertices, faces);

    // Show the gui
    polyscope::state::userCallback = myCallback; // specify the callback
    polyscope::show();
}