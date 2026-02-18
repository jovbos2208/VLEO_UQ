#ifndef SMU_STUBS_H
#define SMU_STUBS_H

#include <Eigen/Dense>

namespace smu {

// Placeholder for smu.rotateAroundPoint
Eigen::Matrix3Xd rotateAroundPoint(const Eigen::Matrix3Xd& points, double angle_rad, const Eigen::Vector3d& direction, const Eigen::Vector3d& point);

// Placeholder for smu.rotateAroundOrigin
Eigen::Matrix3Xd rotateAroundOrigin(const Eigen::Matrix3Xd& points, double angle_rad, const Eigen::Vector3d& direction);

// Placeholder for smu.unitQuat.att.transformVector
Eigen::Vector3d transformVector(const Eigen::Quaterniond& q, const Eigen::Vector3d& v);

// Placeholder for smu.cpm
Eigen::Matrix3d cpm(const Eigen::Vector3d& v);

}

#endif // SMU_STUBS_H
