#ifndef CIRCUMSPHERE_DLL_H
#define CIRCUMSPHERE_DLL_H

#ifdef _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif

#include <vector>

struct Point3D {
    double x, y, z;
};

extern "C" {
    DLL_EXPORT int computeCircumspheres(
        const Point3D* points, const Point3D* normals, int numPoints,
        Point3D* outCenters, double* outRadii, int maxOutputSize, std::string filename);
}

#endif // CIRCUMSPHERE_DLL_H
