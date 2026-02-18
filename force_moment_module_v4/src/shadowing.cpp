#include "shadowing.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace vleo_aerodynamics_core {
namespace {

constexpr double kRayEps = 1e-9;
constexpr double kMinExtent = 1e-3;
constexpr double kExtentMarginRatio = 0.05;
constexpr double kDepthMarginRatio = 0.05;
constexpr int kMinFrameResolution = 128;
constexpr int kMaxFrameResolution = 1024;
constexpr double kPixelsPerTriangle = 6.0;
constexpr double kBarycentricEps = 1e-12;

Eigen::Vector3d computeAxisU(const Eigen::Vector3d& dir_unit) {
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
    if (norm < kRayEps) {
        throw std::runtime_error("Failed to construct perpendicular axis for direction vector.");
    }
    return axis / norm;
}

int determineFrameResolution(int num_triangles) {
    if (const char* env = std::getenv("FM_SHADOW_RES")) {
        if (env && env[0] != '\0') {
            int override_res = std::atoi(env);
            if (override_res > 0) {
                if (override_res < kMinFrameResolution) override_res = kMinFrameResolution;
                if (override_res > kMaxFrameResolution) override_res = kMaxFrameResolution;
                return override_res;
            }
        }
    }

    double pixels_per_triangle = kPixelsPerTriangle;
    if (const char* env = std::getenv("FM_SHADOW_PIXELS_PER_TRIANGLE")) {
        if (env && env[0] != '\0') {
            double ppt = std::atof(env);
            if (ppt > 0.0) {
                pixels_per_triangle = ppt;
            }
        }
    }

    if (num_triangles <= 0) {
        return kMinFrameResolution;
    }
    double estimate = std::sqrt(static_cast<double>(num_triangles)) * pixels_per_triangle;
    int res = static_cast<int>(std::round(estimate));
    if (res < kMinFrameResolution) res = kMinFrameResolution;
    if (res > kMaxFrameResolution) res = kMaxFrameResolution;
    return res;
}

} // namespace

ShadowingResult determineShadowedTriangles(const Eigen::Matrix3Xd& vertices,
                                           const std::vector<Eigen::Vector3d>& centroids,
                                           const std::vector<Eigen::Vector3d>& normals,
                                           const Eigen::Vector3d& dir,
                                           ShadowPixelFrame* pixel_frame_out) {
    const int num_triangles = vertices.cols() / 3;
    if (vertices.cols() % 3 != 0) {
        throw std::runtime_error("Vertices matrix does not contain a whole number of triangles.");
    }
    if (centroids.size() != static_cast<std::size_t>(num_triangles) ||
        normals.size() != static_cast<std::size_t>(num_triangles)) {
        throw std::runtime_error("Mismatch between triangle metadata and vertices.");
    }

    const double dir_norm = dir.norm();
    if (dir_norm < kRayEps) {
        throw std::runtime_error("Direction vector must be non-zero.");
    }
    const Eigen::Vector3d dir_unit = dir / dir_norm;
    const Eigen::Vector3d axis_u = computeAxisU(dir_unit);
    const Eigen::Vector3d axis_v = dir_unit.cross(axis_u);

    const int num_vertices = vertices.cols();
    std::vector<double> vertex_u(num_vertices);
    std::vector<double> vertex_v(num_vertices);
    std::vector<double> vertex_w(num_vertices);

    double min_u = std::numeric_limits<double>::infinity();
    double max_u = -std::numeric_limits<double>::infinity();
    double min_v = std::numeric_limits<double>::infinity();
    double max_v = -std::numeric_limits<double>::infinity();
    double min_w = std::numeric_limits<double>::infinity();
    double max_w = -std::numeric_limits<double>::infinity();

    for (int idx = 0; idx < num_vertices; ++idx) {
        const Eigen::Vector3d vertex = vertices.col(idx);
        const double u = axis_u.dot(vertex);
        const double v = axis_v.dot(vertex);
        const double w = dir_unit.dot(vertex);

        vertex_u[idx] = u;
        vertex_v[idx] = v;
        vertex_w[idx] = w;

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
    const double u_max_frame = max_u + u_margin;
    const double v_min_frame = min_v - v_margin;
    const double v_max_frame = max_v + v_margin;

    const int frame_resolution = determineFrameResolution(num_triangles);
    const double pixel_size_u = (u_max_frame - u_min_frame) / static_cast<double>(frame_resolution);
    const double pixel_size_v = (v_max_frame - v_min_frame) / static_cast<double>(frame_resolution);

    std::vector<int> pixel_triangles(frame_resolution * frame_resolution, -1);
    std::vector<double> pixel_depth(frame_resolution * frame_resolution,
                                    std::numeric_limits<double>::infinity());
    std::vector<uint8_t> triangle_visible(num_triangles, 0);

    for (int tri_idx = 0; tri_idx < num_triangles; ++tri_idx) {
        const bool front_facing = (dir_unit.dot(normals[tri_idx]) <= 0.0);
        if (!front_facing) {
            continue;
        }

        const int base = 3 * tri_idx;
        const double u0 = vertex_u[base + 0];
        const double u1 = vertex_u[base + 1];
        const double u2 = vertex_u[base + 2];
        const double v0 = vertex_v[base + 0];
        const double v1 = vertex_v[base + 1];
        const double v2 = vertex_v[base + 2];
        const double w0 = vertex_w[base + 0];
        const double w1 = vertex_w[base + 1];
        const double w2 = vertex_w[base + 2];

        const double px0 = (u0 - u_min_frame) / pixel_size_u;
        const double px1 = (u1 - u_min_frame) / pixel_size_u;
        const double px2 = (u2 - u_min_frame) / pixel_size_u;
        const double py0 = (v0 - v_min_frame) / pixel_size_v;
        const double py1 = (v1 - v_min_frame) / pixel_size_v;
        const double py2 = (v2 - v_min_frame) / pixel_size_v;

        const double area = (px1 - px0) * (py2 - py0) - (px2 - px0) * (py1 - py0);
        if (std::abs(area) < kBarycentricEps) {
            continue;
        }
        const double inv_area = 1.0 / area;

        const double min_px_center = std::min({px0, px1, px2});
        const double max_px_center = std::max({px0, px1, px2});
        const double min_py_center = std::min({py0, py1, py2});
        const double max_py_center = std::max({py0, py1, py2});

        int px_start = std::max(0, static_cast<int>(std::floor(min_px_center - 0.5)));
        int px_end = std::min(frame_resolution - 1, static_cast<int>(std::ceil(max_px_center - 0.5)));
        int py_start = std::max(0, static_cast<int>(std::floor(min_py_center - 0.5)));
        int py_end = std::min(frame_resolution - 1, static_cast<int>(std::ceil(max_py_center - 0.5)));

        if (px_start > px_end || py_start > py_end) {
            continue;
        }

        for (int py = py_start; py <= py_end; ++py) {
            const double sample_y = static_cast<double>(py) + 0.5;
            for (int px = px_start; px <= px_end; ++px) {
                const double sample_x = static_cast<double>(px) + 0.5;

                const double bary0 = ((px1 - sample_x) * (py2 - sample_y) -
                                      (px2 - sample_x) * (py1 - sample_y)) * inv_area;
                const double bary1 = ((px2 - sample_x) * (py0 - sample_y) -
                                      (px0 - sample_x) * (py2 - sample_y)) * inv_area;
                const double bary2 = 1.0 - bary0 - bary1;

                if (bary0 < -kBarycentricEps ||
                    bary1 < -kBarycentricEps ||
                    bary2 < -kBarycentricEps) {
                    continue;
                }

                const double depth = bary0 * w0 + bary1 * w1 + bary2 * w2;
                const double distance = depth - plane_w;
                if (distance <= kRayEps) {
                    continue;
                }

                const int pixel_index = py * frame_resolution + px;
                if (distance < pixel_depth[pixel_index]) {
                    pixel_depth[pixel_index] = distance;
                    pixel_triangles[pixel_index] = tri_idx;
                }
            }
        }
    }

    std::vector<int> visible_triangles;
    visible_triangles.reserve(num_triangles);
    std::vector<bool> ind_shadowed(num_triangles, true);
    std::vector<int> first_occluder(num_triangles, -1);
    for (int idx = 0; idx < static_cast<int>(pixel_triangles.size()); ++idx) {
        const int triangle_id = pixel_triangles[idx];
        if (triangle_id >= 0 && !triangle_visible[triangle_id]) {
            triangle_visible[triangle_id] = 1;
            visible_triangles.push_back(triangle_id);
            ind_shadowed[triangle_id] = false;
        }
    }

    if (pixel_frame_out) {
        pixel_frame_out->resolution = frame_resolution;
        pixel_frame_out->u_min = u_min_frame;
        pixel_frame_out->v_min = v_min_frame;
        pixel_frame_out->pixel_size_u = pixel_size_u;
        pixel_frame_out->pixel_size_v = pixel_size_v;
        pixel_frame_out->plane_w = plane_w;
        pixel_frame_out->axis_u = axis_u;
        pixel_frame_out->axis_v = axis_v;
        pixel_frame_out->dir_unit = dir_unit;
        pixel_frame_out->pixel_triangle_ids = pixel_triangles;
    }

    return ShadowingResult{std::move(ind_shadowed),
                           std::move(first_occluder),
                           std::move(visible_triangles)};
}

}
