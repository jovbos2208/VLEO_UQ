#include "smu_stubs.h"

namespace smu {

// THIS IS A PLACEHOLDER IMPLEMENTATION
// The user needs to provide the actual implementation of these functions

Eigen::Matrix3Xd rotateAroundPoint(const Eigen::Matrix3Xd& points, double angle_rad, const Eigen::Vector3d& direction, const Eigen::Vector3d& point) {
    Eigen::Vector3d dir = direction;
    if (dir.norm() == 0.0) {
        dir = Eigen::Vector3d::UnitX();
    }
    Eigen::AngleAxisd aa(angle_rad, dir.normalized());
    Eigen::Matrix3d R = aa.toRotationMatrix();

    Eigen::Matrix3Xd out(3, points.cols());
    for (Eigen::Index i = 0; i < points.cols(); ++i) {
        Eigen::Vector3d p = points.col(i) - point;
        out.col(i) = R * p + point;
    }
    return out;
}

Eigen::Matrix3Xd rotateAroundOrigin(const Eigen::Matrix3Xd& points, double angle_rad, const Eigen::Vector3d& direction) {
    Eigen::Vector3d dir = direction;
    if (dir.norm() == 0.0) {
        dir = Eigen::Vector3d::UnitX();
    }
    Eigen::AngleAxisd aa(angle_rad, dir.normalized());
    Eigen::Matrix3d R = aa.toRotationMatrix();

    Eigen::Matrix3Xd out(3, points.cols());
    for (Eigen::Index i = 0; i < points.cols(); ++i) {
        out.col(i) = R * points.col(i);
    }
    return out;
}

Eigen::Vector3d transformVector(const Eigen::Quaterniond& q, const Eigen::Vector3d& v) {
    return q.toRotationMatrix() * v;
}

Eigen::Matrix3d cpm(const Eigen::Vector3d& v) {
    Eigen::Matrix3d m;
    m << 0, -v.z(), v.y(),
         v.z(), 0, -v.x(),
         -v.y(), v.x(), 0;
    return m;
}

}
