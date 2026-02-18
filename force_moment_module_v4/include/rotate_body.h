#ifndef ROTATE_BODY_H
#define ROTATE_BODY_H

#include <Eigen/Dense>
#include <vector>

namespace vleo_aerodynamics_core {

void rotateBody(Eigen::Matrix3Xd& vertices,
                std::vector<Eigen::Vector3d>& centroids,
                std::vector<Eigen::Vector3d>& normals,
                double rotation_angle_rad,
                const Eigen::Vector3d& rotation_direction,
                const Eigen::Vector3d& rotation_hinge_point);

}

#endif // ROTATE_BODY_H
