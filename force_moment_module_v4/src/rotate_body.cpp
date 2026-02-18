#include "rotate_body.h"
#include "smu_stubs.h"

namespace vleo_aerodynamics_core {

void rotateBody(Eigen::Matrix3Xd& vertices,
                std::vector<Eigen::Vector3d>& centroids,
                std::vector<Eigen::Vector3d>& normals,
                double rotation_angle_rad,
                const Eigen::Vector3d& rotation_direction,
                const Eigen::Vector3d& rotation_hinge_point) {

    // Rotate vertices
    vertices = smu::rotateAroundPoint(vertices, rotation_angle_rad, rotation_direction, rotation_hinge_point);

    // Rotate centroids
    Eigen::Matrix3Xd centroids_matrix(3, centroids.size());
    for (size_t i = 0; i < centroids.size(); ++i) {
        centroids_matrix.col(i) = centroids[i];
    }
    centroids_matrix = smu::rotateAroundPoint(centroids_matrix, rotation_angle_rad, rotation_direction, rotation_hinge_point);
    for (size_t i = 0; i < centroids.size(); ++i) {
        centroids[i] = centroids_matrix.col(i);
    }

    // Rotate normals
    Eigen::Matrix3Xd normals_matrix(3, normals.size());
    for (size_t i = 0; i < normals.size(); ++i) {
        normals_matrix.col(i) = normals[i];
    }
    normals_matrix = smu::rotateAroundOrigin(normals_matrix, rotation_angle_rad, rotation_direction);
    for (size_t i = 0; i < normals.size(); ++i) {
        normals[i] = normals_matrix.col(i);
    }
}

}
