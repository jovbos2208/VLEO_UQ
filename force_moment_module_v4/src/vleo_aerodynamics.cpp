#include "vleo_aerodynamics.h"
#include "smu_stubs.h"
#include "shadowing.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <numeric>
#include <utility>

// --- export stuff ---
#include <fstream>

void save_to_VTK(const std::string& filename,
                 const Eigen::Matrix3Xd& vertices,
                 int total_num_faces,
                 const std::vector<bool>& shadow_mask,
                 const std::vector<Eigen::Vector3d>& pressures){
    std::ofstream vtkFile(filename);
    if (!vtkFile.is_open()) {
        std::cerr << "Error: could not open file " << filename << " for writing.\n";
        return;
    }

    const int numVerts = vertices.cols();

    // --- Header ---
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Triangular mesh with scalar cell forces\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET UNSTRUCTURED_GRID\n";

    // --- Points section ---
    vtkFile << "POINTS " << numVerts << " float\n";
    for (int i = 0; i < numVerts; ++i) {
        const Eigen::Vector3d v = vertices.col(i);
        vtkFile << v.x() << " " << v.y() << " " << v.z() << "\n";
    }

    // --- Cells section ---
    vtkFile << "CELLS " << total_num_faces << " " << (total_num_faces*4) <<"\n";
    for (int i = 0; i + 2 < numVerts; i += 3) {
        int i0 = i + 0;
        int i1 = i + 1;
        int i2 = i + 2;
        vtkFile << "3 " << i0 << " " << i1 << " " << i2 << "\n";
    }

    // --- Cell types section ---
    vtkFile << "CELL_TYPES " << total_num_faces <<"\n";
    for (int i = 0; i<total_num_faces; i++){
        vtkFile << "5\n";
    }

    // --- Cell data section ---
    vtkFile << "CELL_DATA " << total_num_faces << "\n";
    vtkFile << "SCALARS shadowed double 1\n";
    vtkFile << "LOOKUP_TABLE default\n";
    for (int i=0; i<total_num_faces;i++){
        vtkFile << (shadow_mask[i] ? 1 : 0) << "\n";
    }

    vtkFile << "VECTORS pressure double\n";
    for (int i = 0; i < total_num_faces; ++i) {
        const Eigen::Vector3d vec = (i < static_cast<int>(pressures.size())) ? pressures[i] : Eigen::Vector3d::Zero();
        vtkFile << vec.x() << " " << vec.y() << " " << vec.z() << "\n";
    }

    vtkFile.close();
    std::cout << "Saved VTK to " << filename << "\n";
}
// --------------------

namespace {

using vleo_aerodynamics_core::Body;

struct VleoWorkingBuffers {
    Eigen::Matrix3Xd all_vertices;
    std::vector<Eigen::Vector3d> all_centroids;
    std::vector<Eigen::Vector3d> all_normals;
    std::vector<double> all_areas;
    std::vector<double> all_temperatures;
    std::vector<double> all_energy_accommodation_coefficients;

    Eigen::VectorXd areas_visible;
    Eigen::VectorXd surface_temperatures_visible;
    Eigen::VectorXd energy_coeffs_visible;
    Eigen::Matrix3Xd normals_visible;
    Eigen::Matrix3Xd centroids_visible;
    Eigen::Matrix3Xd v_rels_visible;

    std::vector<Eigen::Vector3d> pressures_visible;
    std::vector<Eigen::Vector3d> pressures_all_faces;
    std::vector<bool> shadow_mask;
};

VleoWorkingBuffers& workingBuffers(size_t total_vertex_cols, size_t total_num_faces) {
    thread_local VleoWorkingBuffers buffers;
    if (buffers.all_vertices.cols() != static_cast<int>(total_vertex_cols)) {
        buffers.all_vertices.resize(3, static_cast<int>(total_vertex_cols));
    }
    buffers.all_centroids.resize(total_num_faces);
    buffers.all_normals.resize(total_num_faces);
    buffers.all_areas.resize(total_num_faces);
    buffers.all_temperatures.resize(total_num_faces);
    buffers.all_energy_accommodation_coefficients.resize(total_num_faces);
    return buffers;
}

void writeRotatedBody(const Body& body,
                      double rotation_angle_rad,
                      Eigen::Matrix3Xd& vertices_out,
                      int vertex_col_offset,
                      std::vector<Eigen::Vector3d>& centroids_out,
                      size_t centroid_offset,
                      std::vector<Eigen::Vector3d>& normals_out) {
    const int num_vertices = static_cast<int>(body.vertices_B.cols());
    const size_t num_faces = body.centroids_B.size();

    const double dir_norm = body.rotation_direction_B.norm();
    const bool requires_rotation = (std::abs(rotation_angle_rad) > 1e-12) && (dir_norm > 1e-12);

    if (!requires_rotation) {
        vertices_out.block(0, vertex_col_offset, 3, num_vertices) = body.vertices_B;
        for (size_t idx = 0; idx < num_faces; ++idx) {
            centroids_out[centroid_offset + idx] = body.centroids_B[idx];
            normals_out[centroid_offset + idx] = body.normals_B[idx];
        }
        return;
    }

    Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
    Eigen::Vector3d hinge = body.rotation_hinge_point_B;
    Eigen::Vector3d axis = body.rotation_direction_B;
    if (requires_rotation) {
        axis /= dir_norm;
        R = Eigen::AngleAxisd(rotation_angle_rad, axis).toRotationMatrix();
    }

    for (int col = 0; col < num_vertices; ++col) {
        Eigen::Vector3d transformed = body.vertices_B.col(col);
        if (requires_rotation) {
            transformed = R * (transformed - hinge) + hinge;
        }
        vertices_out.col(vertex_col_offset + col) = transformed;
    }

    for (size_t idx = 0; idx < num_faces; ++idx) {
        Eigen::Vector3d centroid = body.centroids_B[idx];
        Eigen::Vector3d normal = body.normals_B[idx];
        if (requires_rotation) {
            centroid = R * (centroid - hinge) + hinge;
            normal = R * normal;
        }
        centroids_out[centroid_offset + idx] = centroid;
        normals_out[centroid_offset + idx] = normal;
    }
}

} // namespace

namespace vleo_aerodynamics_core {

AeroForceAndTorque vleoAerodynamics(const Eigen::Quaterniond& attitude_quaternion_BI,
                                    const Eigen::Vector3d& rotational_velocity_BI_B,
                                    const Eigen::Vector3d& velocity_I_I,
                                    const Eigen::Vector3d& wind_velocity_I_I,
                                    double density,
                                    double temperature,
                                    double particles_mass,
                                    const std::vector<Body>& bodies,
                                    const std::vector<double>& bodies_rotation_angles_rad,
                                    int temperature_ratio_method,
                                    const std::string* vtk_output_path) {

    Eigen::Vector3d v_rel_I = wind_velocity_I_I - velocity_I_I;
    Eigen::Vector3d v_rel_B = smu::transformVector(attitude_quaternion_BI, v_rel_I);
    
    std::string debug_vtk_path;
    if (const char* env = std::getenv("FM_DEBUG_VTK"); env && env[0] != '\0') {
        debug_vtk_path = env;
    }
    const bool debug_vtk_enabled = false;                                    

    size_t total_num_faces = 0;
    for (const auto& body : bodies) {
        total_num_faces += body.vertices_B.cols() / 3;
    }

    // Preallocate vertices storage to avoid repeated conservativeResize
    size_t total_vertex_cols = 0;
    for (const auto& body : bodies) {
        total_vertex_cols += static_cast<size_t>(body.vertices_B.cols());
    }
    VleoWorkingBuffers& buffers = workingBuffers(total_vertex_cols, total_num_faces);
    Eigen::Matrix3Xd& all_vertices = buffers.all_vertices;
    auto& all_centroids = buffers.all_centroids;
    auto& all_normals = buffers.all_normals;
    auto& all_areas = buffers.all_areas;
    auto& all_temperatures = buffers.all_temperatures;
    auto& all_energy_accommodation_coefficients = buffers.all_energy_accommodation_coefficients;
    
    size_t face_idx = 0;
    int vertex_col_offset = 0;
    for (size_t i = 0; i < bodies.size(); ++i) {
        const Body& body = bodies[i];
        const double rotation_angle = (i < bodies_rotation_angles_rad.size()) ? bodies_rotation_angles_rad[i] : 0.0;
        writeRotatedBody(body,
                         rotation_angle,
                         all_vertices,
                         vertex_col_offset,
                         all_centroids,
                         face_idx,
                         all_normals);

        vertex_col_offset += static_cast<int>(body.vertices_B.cols());

        for (size_t j = 0; j < body.centroids_B.size(); ++j) {
            all_areas[face_idx] = body.areas[j];
            all_temperatures[face_idx] = body.temperatures_K[j];
            all_energy_accommodation_coefficients[face_idx] = body.energy_accommodation_coefficients[j];
            face_idx++;
        }
    }

    Eigen::Vector3d v_rel_dir_B = v_rel_B.normalized();
    std::vector<int> visible_triangles;
    if (const char* env = std::getenv("FM_DISABLE_SHADOWING"); env && env[0] != '\0' && std::string(env) != "0") {
        visible_triangles.resize(total_num_faces);
        std::iota(visible_triangles.begin(), visible_triangles.end(), 0);
    } else {
        ShadowPixelFrame pixel_frame;
        ShadowingResult shadowing = determineShadowedTriangles(all_vertices,
                                                               all_centroids,
                                                               all_normals,
                                                               v_rel_dir_B,
                                                               debug_vtk_enabled ? &pixel_frame : nullptr);
        visible_triangles = std::move(shadowing.visible_triangles);
    }

    // Gather visible faces into contiguous buffers
    const size_t visible_count = visible_triangles.size();
    buffers.areas_visible.resize(visible_count);
    buffers.surface_temperatures_visible.resize(visible_count);
    buffers.energy_coeffs_visible.resize(visible_count);
    buffers.normals_visible.resize(3, visible_count);
    buffers.centroids_visible.resize(3, visible_count);
    buffers.v_rels_visible.resize(3, visible_count);

    for (size_t idx = 0; idx < visible_count; ++idx) {
        const int tri_idx = visible_triangles[idx];
        buffers.areas_visible(idx) = all_areas[tri_idx];
        buffers.surface_temperatures_visible(idx) = all_temperatures[tri_idx];
        buffers.energy_coeffs_visible(idx) = all_energy_accommodation_coefficients[tri_idx];
        buffers.normals_visible.col(idx) = all_normals[tri_idx];
        buffers.centroids_visible.col(idx) = all_centroids[tri_idx];
        buffers.v_rels_visible.col(idx) = v_rel_B - rotational_velocity_BI_B.cross(all_centroids[tri_idx]);
    }

    auto& pressures_visible = buffers.pressures_visible;
    auto& pressures_all_faces = buffers.pressures_all_faces;
    auto& shadow_mask = buffers.shadow_mask;
    if (vtk_output_path) {
        pressures_visible.clear();
        pressures_visible.reserve(visible_count);
        pressures_all_faces.assign(total_num_faces, Eigen::Vector3d::Zero());
        shadow_mask.assign(total_num_faces, true);
        for (int tri_idx : visible_triangles) {
            shadow_mask[tri_idx] = false;
        }
    }

    AeroForceAndTorque result = calcAeroForceAndTorque(buffers.areas_visible,
                                                       buffers.normals_visible,
                                                       buffers.centroids_visible,
                                                       buffers.v_rels_visible,
                                                       density,
                                                       temperature,
                                                       buffers.surface_temperatures_visible,
                                                       buffers.energy_coeffs_visible,
                                                       particles_mass,
                                                       temperature_ratio_method,
                                                       vtk_output_path ? &pressures_visible : nullptr);

    if (vtk_output_path) {
        for (size_t idx = 0; idx < visible_triangles.size(); ++idx) {
            pressures_all_faces[visible_triangles[idx]] = pressures_visible[idx];
        }
        save_to_VTK(*vtk_output_path,
                    all_vertices,
                    static_cast<int>(total_num_faces),
                    shadow_mask,
                    pressures_all_faces);
    }

    return result;
}

}
