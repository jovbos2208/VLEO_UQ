#include "body_importer.h"
#include "vleo_aerodynamics.h"

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

using json = nlohmann::json;

int main() {
    std::string input;
    if (!std::getline(std::cin, input)) {
        std::cerr << "No input provided on stdin.\n";
        return 1;
    }

    json payload;
    try {
        payload = json::parse(input);
    } catch (const std::exception& e) {
        std::cerr << "Failed to parse JSON input: " << e.what() << "\n";
        return 1;
    }

    if (!payload.contains("obj_path") || !payload.contains("data")) {
        std::cerr << "JSON must contain 'obj_path' and 'data' fields.\n";
        return 1;
    }

    const std::string obj_path = payload.at("obj_path").get<std::string>();
    if (!std::filesystem::exists(obj_path)) {
        std::cerr << "OBJ file not found: " << obj_path << "\n";
        return 1;
    }

    const std::vector<double> data = payload.at("data").get<std::vector<double>>();
    constexpr int kStride = 6; // Vx Vy Vz rho T_K s
    if (data.empty() || (data.size() % kStride != 0)) {
        std::cerr << "Input data must be a flat list with length divisible by "
                  << kStride << ".\n";
        return 1;
    }

    Eigen::Map<const Eigen::Matrix<double, kStride, Eigen::Dynamic>> input_matrix(
        data.data(), kStride, data.size() / kStride);

    using namespace vleo_aerodynamics_core;

    std::vector<std::string> object_files = {obj_path};
    Eigen::Matrix3Xd hinge_points(3, 1);
    Eigen::Matrix3Xd hinge_axes(3, 1);
    hinge_points.setZero();
    hinge_axes.col(0) = Eigen::Vector3d::UnitZ();
    std::vector<std::vector<double>> temperatures_K = {std::vector<double>{300.0}};
    std::vector<std::vector<double>> eac = {std::vector<double>{0.9}};

    Eigen::Matrix3d dcm = Eigen::Matrix3d::Identity();
    Eigen::Vector3d com = Eigen::Vector3d::Zero();
    const std::vector<Body> bodies = importMultipleBodies(
        object_files, hinge_points, hinge_axes, temperatures_K, eac, dcm, com);
    const std::vector<double> angles = {0.0};

    const Eigen::Quaterniond q_BI(1.0, 0.0, 0.0, 0.0);
    const Eigen::Vector3d omega_BI_B = Eigen::Vector3d::Zero();
    const Eigen::Vector3d wind_I = Eigen::Vector3d::Zero();
    const double kB = 1.38064852e-23;
    const int temp_ratio_method = 2;

    for (int c = 0; c < input_matrix.cols(); ++c) {
        const double Vx = input_matrix(0, c);
        const double Vy = input_matrix(1, c);
        const double Vz = input_matrix(2, c);
        const double rho = input_matrix(3, c);
        const double T_K = input_matrix(4, c);
        const double s = input_matrix(5, c);

        Eigen::Vector3d force = Eigen::Vector3d::Zero();
        Eigen::Vector3d torque = Eigen::Vector3d::Zero();

        const Eigen::Vector3d V_I(Vx, Vy, Vz);
        const double V_mag = V_I.norm();
        if (V_mag > 0.0 && s > 0.0) {
            const double ratio = s / V_mag;
            const double m_particle = 2.0 * kB * T_K * ratio * ratio;
            const auto out = vleoAerodynamics(
                q_BI, omega_BI_B, V_I, wind_I, rho, T_K, m_particle, bodies, angles,
                temp_ratio_method, nullptr);
            force = out.force;
            torque = out.torque;
        }

        std::cout << force.transpose() << " " << torque.transpose() << '\n';
    }
    return 0;
}
