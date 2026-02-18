#ifndef BODY_IMPORTER_H
#define BODY_IMPORTER_H

#include "obj_loader.h"
#include <Eigen/Dense>
#include <string>
#include <vector>

namespace vleo_aerodynamics_core {

struct Body {
    Eigen::Matrix3Xd vertices_B;
    std::vector<Eigen::Vector3d> centroids_B;
    std::vector<Eigen::Vector3d> normals_B;
    std::vector<double> areas;
    Eigen::Vector3d rotation_hinge_point_B;
    Eigen::Vector3d rotation_direction_B;
    std::vector<double> temperatures_K;
    std::vector<double> energy_accommodation_coefficients;
};

std::vector<Body> importMultipleBodies(const std::vector<std::string>& object_files,
                                       const Eigen::Matrix3Xd& rotation_hinge_points_CAD,
                                       const Eigen::Matrix3Xd& rotation_directions_CAD,
                                       const std::vector<std::vector<double>>& temperatures_K,
                                       const std::vector<std::vector<double>>& energy_accommodation_coefficients,
                                       const Eigen::Matrix3d& DCM_B_from_CAD,
                                       const Eigen::Vector3d& CoM_CAD);

}

#endif // BODY_IMPORTER_H
