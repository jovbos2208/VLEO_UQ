#ifndef SHADOWING_H
#define SHADOWING_H

#include <Eigen/Dense>
#include <vector>

namespace vleo_aerodynamics_core {

struct ShadowingResult {
    std::vector<bool> is_shadowed;
    std::vector<int> first_occluder;
    std::vector<int> visible_triangles;
};

struct ShadowPixelFrame {
    int resolution = 0;
    double u_min = 0.0;
    double v_min = 0.0;
    double pixel_size_u = 0.0;
    double pixel_size_v = 0.0;
    double plane_w = 0.0;
    Eigen::Vector3d axis_u = Eigen::Vector3d::Zero();
    Eigen::Vector3d axis_v = Eigen::Vector3d::Zero();
    Eigen::Vector3d dir_unit = Eigen::Vector3d::Zero();
    std::vector<int> pixel_triangle_ids;
};

ShadowingResult determineShadowedTriangles(const Eigen::Matrix3Xd& vertices,
                                           const std::vector<Eigen::Vector3d>& centroids,
                                           const std::vector<Eigen::Vector3d>& normals,
                                           const Eigen::Vector3d& dir,
                                           ShadowPixelFrame* pixel_frame_out = nullptr);

}

#endif // SHADOWING_H
