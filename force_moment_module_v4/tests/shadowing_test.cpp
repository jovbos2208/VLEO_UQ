#include "shadowing.h"
#include "body_importer.h"
#include "rotate_body.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#ifndef FORCE_MOMENT_MODULE_SOURCE_DIR
#define FORCE_MOMENT_MODULE_SOURCE_DIR "."
#endif

namespace vleo_aerodynamics_core {

struct MeshBuffers {
    Eigen::Matrix3Xd vertices;
    std::vector<Eigen::Vector3d> centroids;
    std::vector<Eigen::Vector3d> normals;
};

static MeshBuffers makeStackedTriangles() {
    MeshBuffers mesh;
    mesh.vertices.resize(3, 6);
    mesh.centroids.resize(2);
    mesh.normals.resize(2, Eigen::Vector3d::UnitX());

    Eigen::Matrix<double, 3, 3> front_triangle;
    front_triangle.col(0) = Eigen::Vector3d(0.0, 0.0, 0.0);
    front_triangle.col(1) = Eigen::Vector3d(0.0, 1.0, 0.0);
    front_triangle.col(2) = Eigen::Vector3d(0.0, 0.0, 1.0);

    Eigen::Matrix<double, 3, 3> rear_triangle;
    rear_triangle.col(0) = Eigen::Vector3d(-1.0, 0.0, 0.0);
    rear_triangle.col(1) = Eigen::Vector3d(-1.0, 1.0, 0.0);
    rear_triangle.col(2) = Eigen::Vector3d(-1.0, 0.0, 1.0);

    mesh.vertices.block(0, 0, 3, 3) = front_triangle;
    mesh.vertices.block(0, 3, 3, 3) = rear_triangle;

    mesh.centroids[0] = front_triangle.rowwise().mean();
    mesh.centroids[1] = rear_triangle.rowwise().mean();
    return mesh;
}

static Eigen::Matrix<double, 3, 3> extractTriangle(const Eigen::Matrix3Xd& vertices, int idx) {
    return vertices.block(0, 3 * idx, 3, 3);
}

static Eigen::Vector3d computeAxisU(const Eigen::Vector3d& dir_unit) {
    Eigen::Vector3d helper;
    if (std::abs(dir_unit.x()) <= std::abs(dir_unit.y()) &&
        std::abs(dir_unit.x()) <= std::abs(dir_unit.z())) {
        helper = Eigen::Vector3d::UnitX();
    } else if (std::abs(dir_unit.y()) <= std::abs(dir_unit.x()) &&
               std::abs(dir_unit.y()) <= std::abs(dir_unit.z())) {
        helper = Eigen::Vector3d::UnitY();
    } else {
        helper = Eigen::Vector3d::UnitZ();
    }
    Eigen::Vector3d axis = helper - dir_unit * dir_unit.dot(helper);
    double norm = axis.norm();
    if (norm < 1e-9) {
        throw std::runtime_error("Failed to build perpendicular axis for reference shadowing.");
    }
    return axis / norm;
}

static bool rayIntersectsTriangle(const Eigen::Vector3d& origin,
                                  const Eigen::Vector3d& direction,
                                  const Eigen::Matrix<double, 3, 3>& triangle,
                                  double& out_distance) {
    constexpr double EPS = 1e-9;
    const Eigen::Vector3d edge1 = triangle.col(1) - triangle.col(0);
    const Eigen::Vector3d edge2 = triangle.col(2) - triangle.col(0);
    const Eigen::Vector3d h = direction.cross(edge2);
    const double a = edge1.dot(h);
    if (std::abs(a) < EPS) {
        return false;
    }
    const double f = 1.0 / a;
    const Eigen::Vector3d s = origin - triangle.col(0);
    const double u = f * s.dot(h);
    if (u < 0.0 || u > 1.0) {
        return false;
    }
    const Eigen::Vector3d q = s.cross(edge1);
    const double v = f * direction.dot(q);
    if (v < 0.0 || (u + v) > 1.0) {
        return false;
    }
    const double t = f * edge2.dot(q);
    if (t <= EPS) {
        return false;
    }
    out_distance = t;
    return true;
}

static std::vector<bool> referenceShadowing(const MeshBuffers& mesh,
                                            const Eigen::Vector3d& dir) {
    constexpr double kMinExtent = 1e-3;
    constexpr double kExtentMarginRatio = 0.05;
    constexpr double kDepthMarginRatio = 0.05;
    constexpr int kResolution = 512;

    const Eigen::Vector3d dir_unit = dir.normalized();
    const Eigen::Vector3d axis_u = computeAxisU(dir_unit);
    const Eigen::Vector3d axis_v = dir_unit.cross(axis_u);

    double min_u = std::numeric_limits<double>::infinity();
    double max_u = -std::numeric_limits<double>::infinity();
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();
    double min_w = std::numeric_limits<double>::infinity();
    double max_w = -std::numeric_limits<double>::infinity();

    for (int idx = 0; idx < mesh.vertices.cols(); ++idx) {
        const Eigen::Vector3d vertex = mesh.vertices.col(idx);
        const double u = axis_u.dot(vertex);
        const double v = axis_v.dot(vertex);
        const double w = dir_unit.dot(vertex);
        min_u = std::min(min_u, u);
        max_u = std::max(max_u, u);
        min_v = std::min(min_v, v);
        max_v = std::max(max_v, v);
        min_w = std::min(min_w, w);
        max_w = std::max(max_w, w);
    }

    const double depth_range = std::max(max_w - min_w, kMinExtent);
    const double extent_u = std::max(max_u - min_u, kMinExtent);
    const double extent_v = std::max(max_v - min_v, kMinExtent);
    const double u_margin = std::max(kExtentMarginRatio * extent_u, kMinExtent);
    const double v_margin = std::max(kExtentMarginRatio * extent_v, kMinExtent);
    const double plane_w = min_w - std::max(kDepthMarginRatio * depth_range, kMinExtent);

    const double u_min_frame = min_u - u_margin;
    const double v_min_frame = min_v - v_margin;
    const double pixel_size_u = (extent_u + 2.0 * u_margin) / static_cast<double>(kResolution);
    const double pixel_size_v = (extent_v + 2.0 * v_margin) / static_cast<double>(kResolution);

    const std::size_t num_triangles = mesh.centroids.size();
    std::vector<bool> triangle_hit(num_triangles, false);
    const Eigen::Vector3d ray_dir = dir_unit;

    std::vector<Eigen::Matrix<double, 3, 3>> triangles(num_triangles);
    for (std::size_t i = 0; i < num_triangles; ++i) {
        triangles[i] = extractTriangle(mesh.vertices, static_cast<int>(i));
    }

    for (int py = 0; py < kResolution; ++py) {
        const double v_coord = v_min_frame + (py + 0.5) * pixel_size_v;
        for (int px = 0; px < kResolution; ++px) {
            const double u_coord = u_min_frame + (px + 0.5) * pixel_size_u;
            Eigen::Vector3d origin = axis_u * u_coord + axis_v * v_coord + dir_unit * plane_w;

            double best_distance = std::numeric_limits<double>::infinity();
            int best_triangle = -1;
            for (std::size_t t_idx = 0; t_idx < num_triangles; ++t_idx) {
                double hit_distance = 0.0;
                if (rayIntersectsTriangle(origin, ray_dir, triangles[t_idx], hit_distance)) {
                    if (hit_distance < best_distance) {
                        best_distance = hit_distance;
                        best_triangle = static_cast<int>(t_idx);
                    }
                }
            }
            if (best_triangle >= 0) {
                triangle_hit[best_triangle] = true;
            }
        }
    }

    std::vector<bool> result(num_triangles, false);
    for (std::size_t i = 0; i < num_triangles; ++i) {
        const bool front_facing = (dir_unit.dot(mesh.normals[i]) <= 0.0);
        if (!front_facing || !triangle_hit[i]) {
            result[i] = true;
        }
    }
    return result;
}

static MeshBuffers flattenBodies(const std::vector<Body>& bodies,
                                 const std::vector<double>& angles_rad) {
    assert(bodies.size() == angles_rad.size());

    std::size_t total_faces = 0;
    std::size_t total_vertex_cols = 0;
    for (const auto& body : bodies) {
        total_faces += body.centroids_B.size();
        total_vertex_cols += static_cast<std::size_t>(body.vertices_B.cols());
    }

    MeshBuffers mesh;
    mesh.vertices.resize(3, static_cast<int>(total_vertex_cols));
    mesh.centroids.resize(total_faces);
    mesh.normals.resize(total_faces);

    std::size_t face_idx = 0;
    int vertex_offset = 0;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        Eigen::Matrix3Xd rotated_vertices = bodies[i].vertices_B;
        std::vector<Eigen::Vector3d> rotated_centroids = bodies[i].centroids_B;
        std::vector<Eigen::Vector3d> rotated_normals = bodies[i].normals_B;

        rotateBody(rotated_vertices,
                   rotated_centroids,
                   rotated_normals,
                   angles_rad[i],
                   bodies[i].rotation_direction_B,
                   bodies[i].rotation_hinge_point_B);

        mesh.vertices.block(0, vertex_offset, 3, rotated_vertices.cols()) = rotated_vertices;
        vertex_offset += rotated_vertices.cols();

        for (std::size_t j = 0; j < rotated_centroids.size(); ++j) {
            mesh.centroids[face_idx] = rotated_centroids[j];
            mesh.normals[face_idx] = rotated_normals[j];
            ++face_idx;
        }
    }

    return mesh;
}

static MeshBuffers makeShuttlecockMesh() {
    namespace fs = std::filesystem;
    const fs::path root = fs::path(FORCE_MOMENT_MODULE_SOURCE_DIR);
    const fs::path data_dir = root / "data";

    std::vector<std::string> object_files = {
        (data_dir / "shuttle_box.obj").string(),
        (data_dir / "wing_pos_x.obj").string(),
        (data_dir / "wing_neg_x.obj").string(),
        (data_dir / "wing_pos_y.obj").string(),
        (data_dir / "wing_neg_y.obj").string()
    };

    const double box_half = 0.05;
    Eigen::Matrix3Xd hinge(3, object_files.size());
    Eigen::Matrix3Xd axes(3, object_files.size());
    hinge.setZero();
    axes.setZero();

    hinge.col(1) << +box_half, 0.0, 0.0;
    hinge.col(2) << -box_half, 0.0, 0.0;
    hinge.col(3) << 0.0, +box_half, 0.0;
    hinge.col(4) << 0.0, -box_half, 0.0;

    axes.col(0) = Eigen::Vector3d::UnitZ();
    axes.col(1) = Eigen::Vector3d::UnitY();
    axes.col(2) = Eigen::Vector3d::UnitY();
    axes.col(3) = Eigen::Vector3d::UnitX();
    axes.col(4) = Eigen::Vector3d::UnitX();

    std::vector<std::vector<double>> temps_K(object_files.size(), std::vector<double>{300.0});
    std::vector<std::vector<double>> eac(object_files.size(), std::vector<double>{0.9});

    Eigen::Matrix3d DCM = Eigen::Matrix3d::Identity();
    Eigen::Vector3d CoM = Eigen::Vector3d::Zero();

    auto bodies = importMultipleBodies(object_files,
                                       hinge,
                                       axes,
                                       temps_K,
                                       eac,
                                       DCM,
                                       CoM);

    std::vector<double> angles_rad(bodies.size(), 0.0);
    return flattenBodies(bodies, angles_rad);
}

static std::vector<int> visibleTrianglesFromMask(const std::vector<bool>& shadow_mask) {
    std::vector<int> indices;
    indices.reserve(shadow_mask.size());
    for (std::size_t i = 0; i < shadow_mask.size(); ++i) {
        if (!shadow_mask[i]) {
            indices.push_back(static_cast<int>(i));
        }
    }
    return indices;
}

static bool runSimpleOcclusionTest() {
    const MeshBuffers mesh = makeStackedTriangles();
    const Eigen::Vector3d dir = Eigen::Vector3d(-1.0, 0.0, 0.0);

    const auto expected = referenceShadowing(mesh, dir);
    assert(expected.size() == 2);
    assert(expected[0] == false && expected[1] == true);

    const auto actual = determineShadowedTriangles(mesh.vertices, mesh.centroids, mesh.normals, dir);

    if (actual.is_shadowed != expected) {
        std::cerr << "Simple occlusion test failed.\n";
        return false;
    }

    auto expected_visible = visibleTrianglesFromMask(expected);
    auto actual_visible = actual.visible_triangles;
    std::sort(expected_visible.begin(), expected_visible.end());
    std::sort(actual_visible.begin(), actual_visible.end());
    if (actual_visible != expected_visible) {
        std::cerr << "Simple occlusion test failed (visible triangles mismatch).\n";
        return false;
    }
    return true;
}

static bool runShuttlecockTest() {
    const MeshBuffers mesh = makeShuttlecockMesh();
    const Eigen::Vector3d dir = Eigen::Vector3d(-1.0, 0.2, 0.0).normalized();

    const auto expected = referenceShadowing(mesh, dir);
    const auto actual = determineShadowedTriangles(mesh.vertices, mesh.centroids, mesh.normals, dir);

    if (actual.is_shadowed != expected) {
        std::size_t mismatch_shadow = 0;
        for (std::size_t i = 0; i < expected.size(); ++i) {
            if (expected[i] != actual.is_shadowed[i]) {
                ++mismatch_shadow;
            }
        }
        const double mismatch_ratio =
            static_cast<double>(mismatch_shadow) / static_cast<double>(expected.size());
        if (mismatch_ratio > 0.1) {
            std::cerr << "Shuttlecock test failed. Shadow mismatches: "
                      << mismatch_shadow << " (" << mismatch_ratio * 100.0 << "%)\n";
            return false;
        }
    }

    auto expected_visible = visibleTrianglesFromMask(expected);
    auto actual_visible = actual.visible_triangles;
    std::sort(expected_visible.begin(), expected_visible.end());
    std::sort(actual_visible.begin(), actual_visible.end());
    std::size_t mismatch_visible = 0;
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < expected_visible.size() || j < actual_visible.size()) {
        if (j >= actual_visible.size() ||
            (i < expected_visible.size() && expected_visible[i] < actual_visible[j])) {
            ++mismatch_visible;
            ++i;
        } else if (i >= expected_visible.size() ||
                   actual_visible[j] < expected_visible[i]) {
            ++mismatch_visible;
            ++j;
        } else {
            ++i;
            ++j;
        }
    }
    const double mismatch_ratio_visible =
        expected.empty()
            ? 0.0
            : static_cast<double>(mismatch_visible) / static_cast<double>(expected.size());
    if (mismatch_ratio_visible > 0.1) {
        std::cerr << "Shuttlecock test failed (visible triangles mismatch: "
                  << mismatch_visible << " entries).\n";
        return false;
    }
    return true;
}

} // namespace vleo_aerodynamics_core

int main() {
    using namespace vleo_aerodynamics_core;
    bool ok = true;
    ok &= runSimpleOcclusionTest();
    ok &= runShuttlecockTest();
    if (!ok) {
        return 1;
    }
    std::cout << "All shadowing tests passed.\n";
    return 0;
}
