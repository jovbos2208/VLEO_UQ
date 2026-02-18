#ifndef OBJ_LOADER_H
#define OBJ_LOADER_H

#include <string>
#include <vector>
#include <Eigen/Dense>

namespace vleo_aerodynamics_core {

// Struct to hold the result of loading an OBJ file
struct ObjData {
    Eigen::Matrix3Xd vertices;
    std::vector<Eigen::Vector3d> centroids;
    std::vector<Eigen::Vector3d> normals;
    std::vector<double> areas;
};

// Function to load an OBJ file
ObjData loadObjFile(const std::string& filename);

}

#endif // OBJ_LOADER_H
