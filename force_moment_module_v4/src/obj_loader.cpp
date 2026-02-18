#include "obj_loader.h"
#include "tiny_obj_loader.h"
#include <iostream>

namespace vleo_aerodynamics_core {

ObjData loadObjFile(const std::string& filename) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str())) {
        throw std::runtime_error(warn + err);
    }

    ObjData obj_data;
    obj_data.vertices.resize(3, 0); // Initialize with 3 rows and 0 columns

    auto subdivide_and_append = [&](const Eigen::Vector3d& A,
                                    const Eigen::Vector3d& B,
                                    const Eigen::Vector3d& C,
                                    int k) {
        // Uniformly subdivide triangle ABC by splitting each edge into k segments.
        // Generates k^2 smaller triangles. k>=1, k=1 means original triangle.
        if (k <= 1) {
            Eigen::Matrix3d face_vertices_fixed;
            face_vertices_fixed.col(0) = A;
            face_vertices_fixed.col(1) = B;
            face_vertices_fixed.col(2) = C;
            obj_data.vertices.conservativeResize(Eigen::NoChange, obj_data.vertices.cols() + 3);
            obj_data.vertices.rightCols<3>() = face_vertices_fixed;

            Eigen::Vector3d centroid = face_vertices_fixed.rowwise().mean();
            obj_data.centroids.push_back(centroid);
            Eigen::Vector3d v0 = face_vertices_fixed.col(1) - face_vertices_fixed.col(0);
            Eigen::Vector3d v1 = face_vertices_fixed.col(2) - face_vertices_fixed.col(0);
            Eigen::Vector3d normal = v0.cross(v1);
            double area = 0.5 * normal.norm();
            if (area > 0) normal.normalize();
            obj_data.areas.push_back(area);
            obj_data.normals.push_back(normal);
            return;
        }

        // Precompute barycentric grid points P(i,j) where i+j <= k
        std::vector<Eigen::Vector3d> P;
        P.reserve((k + 1) * (k + 2) / 2);
        auto idx = [&](int i, int j) {
            // map (i,j) with i>=0, j>=0, i+j<=k to linear index
            // rows of constant i: j=0..(k-i)
            int offset = 0;
            for (int ii = 0; ii < i; ++ii) offset += (k - ii + 1);
            return offset + j;
        };
        for (int i = 0; i <= k; ++i) {
            for (int j = 0; j <= k - i; ++j) {
                double a = double(k - i - j) / k;
                double b = double(i) / k;
                double c = double(j) / k;
                P.push_back(a * A + b * B + c * C);
            }
        }

        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k - i; ++j) {
                // Small triangle 1: (i,j), (i+1,j), (i,j+1)
                Eigen::Vector3d p0 = P[idx(i, j)];
                Eigen::Vector3d p1 = P[idx(i + 1, j)];
                Eigen::Vector3d p2 = P[idx(i, j + 1)];
                {
                    Eigen::Matrix3d tri;
                    tri.col(0) = p0; tri.col(1) = p1; tri.col(2) = p2;
                    obj_data.vertices.conservativeResize(Eigen::NoChange, obj_data.vertices.cols() + 3);
                    obj_data.vertices.rightCols<3>() = tri;
                    Eigen::Vector3d centroid = tri.rowwise().mean();
                    obj_data.centroids.push_back(centroid);
                    Eigen::Vector3d v0 = tri.col(1) - tri.col(0);
                    Eigen::Vector3d v1 = tri.col(2) - tri.col(0);
                    Eigen::Vector3d normal = v0.cross(v1);
                    double area = 0.5 * normal.norm();
                    if (area > 0) normal.normalize();
                    obj_data.areas.push_back(area);
                    obj_data.normals.push_back(normal);
                }

                // Small triangle 2: (i+1,j), (i+1,j+1), (i,j+1)
                if (i + j < k - 1) {
                    Eigen::Vector3d q0 = P[idx(i + 1, j)];
                    Eigen::Vector3d q1 = P[idx(i + 1, j + 1)];
                    Eigen::Vector3d q2 = P[idx(i, j + 1)];
                    Eigen::Matrix3d tri;
                    tri.col(0) = q0; tri.col(1) = q1; tri.col(2) = q2;
                    obj_data.vertices.conservativeResize(Eigen::NoChange, obj_data.vertices.cols() + 3);
                    obj_data.vertices.rightCols<3>() = tri;
                    Eigen::Vector3d centroid = tri.rowwise().mean();
                    obj_data.centroids.push_back(centroid);
                    Eigen::Vector3d v0 = tri.col(1) - tri.col(0);
                    Eigen::Vector3d v1 = tri.col(2) - tri.col(0);
                    Eigen::Vector3d normal = v0.cross(v1);
                    double area = 0.5 * normal.norm();
                    if (area > 0) normal.normalize();
                    obj_data.areas.push_back(area);
                    obj_data.normals.push_back(normal);
                }
            }
        }
    };

    for (const auto& shape : shapes) {
        for (size_t f = 0; f < shape.mesh.num_face_vertices.size(); ++f) {
            if (shape.mesh.num_face_vertices[f] != 3) {
                // Only triangles are supported
                continue;
            }

            Eigen::Vector3d A, B, C;
            for (size_t v = 0; v < 3; ++v) {
                tinyobj::index_t idxv = shape.mesh.indices[3 * f + v];
                double x = attrib.vertices[3 * idxv.vertex_index + 0];
                double y = attrib.vertices[3 * idxv.vertex_index + 1];
                double z = attrib.vertices[3 * idxv.vertex_index + 2];
                if (v == 0) A = Eigen::Vector3d(x, y, z);
                else if (v == 1) B = Eigen::Vector3d(x, y, z);
                else C = Eigen::Vector3d(x, y, z);
            }
            // Subdivide each input triangle into k^2 smaller triangles. k=3 => ~9x faces (target ~10x)
            const int k = 1;
            subdivide_and_append(A, B, C, k);
        }
    }

    return obj_data;
}

}
