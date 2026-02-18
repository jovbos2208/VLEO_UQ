#include "shuttlecock_aero.h"
#include "vleo_aerodynamics.h"
#include <vector>
#include <string>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <cstdlib>

namespace vleo_aerodynamics_core {

static std::string resolve_path_ancestor(const std::string& name) {
    namespace fs = std::filesystem;
    fs::path here = fs::path(__FILE__).parent_path();
    std::vector<fs::path> candidates = {
        here / ".." / "data" / name,
        fs::current_path() / "data" / name,
        here / name
    };
    for (const auto& cand : candidates) {
        fs::path abs = fs::absolute(cand);
        if (fs::exists(abs)) {
            return abs.string();
        }
    }
    return fs::absolute(here / ".." / "data" / name).string();
}

// should be obsolete
// float positiveValue(float x) {
//     if (x < 0.0f)
//         return 0.0f;
//     return x;
// }

namespace {

struct ShuttlecockGeometryConfig {
    std::vector<std::string> object_files;
    Eigen::Matrix3Xd hinge_points;
    Eigen::Matrix3Xd rotation_axes;
    std::vector<std::vector<double>> temperatures_K;
    std::vector<std::vector<double>> energy_coeffs;
    Eigen::Matrix3d DCM;
    Eigen::Vector3d CoM;
};

const ShuttlecockGeometryConfig& shuttlecockGeometryConfig() {
    static ShuttlecockGeometryConfig config = [] {
        ShuttlecockGeometryConfig cfg;
        cfg.object_files = {
            resolve_path_ancestor("mainBody.obj"),
            resolve_path_ancestor("wing_1.obj"),
            resolve_path_ancestor("wing_2.obj"),
            resolve_path_ancestor("wing_3.obj"),
            resolve_path_ancestor("wing_4.obj"),
        };

        const double box_half = 0.05;
        cfg.hinge_points.resize(3, 5);
        cfg.rotation_axes.resize(3, 5);
        cfg.hinge_points.col(0) << 0.0, 0.0, 0.0;
        cfg.hinge_points.col(1) << -0.34965, 0.0, -box_half;
        cfg.hinge_points.col(2) << -0.34965, -box_half, 0.0;
        cfg.hinge_points.col(3) << -0.34965, 0.0, +box_half;
        cfg.hinge_points.col(4) << -0.34965, +box_half, 0.0;
        cfg.rotation_axes.setZero();
        cfg.rotation_axes.col(0) = Eigen::Vector3d::UnitZ();
        cfg.rotation_axes.col(1) = Eigen::Vector3d::UnitY();
        cfg.rotation_axes.col(2) = Eigen::Vector3d::UnitZ();
        cfg.rotation_axes.col(3) = Eigen::Vector3d::UnitY();
        cfg.rotation_axes.col(4) = Eigen::Vector3d::UnitZ();

        cfg.temperatures_K.assign(5, std::vector<double>{300.0});
        cfg.energy_coeffs.assign(5, std::vector<double>{0.9});
        cfg.DCM = Eigen::Matrix3d::Identity();
        cfg.CoM = Eigen::Vector3d::Zero();
        return cfg;
    }();
    return config;
}

const std::vector<Body>& shuttlecockBaseBodies() {
    static std::vector<Body> base_bodies = [] {
        const auto& cfg = shuttlecockGeometryConfig();
        return importMultipleBodies(cfg.object_files,
                                    cfg.hinge_points,
                                    cfg.rotation_axes,
                                    cfg.temperatures_K,
                                    cfg.energy_coeffs,
                                    cfg.DCM,
                                    cfg.CoM);
    }();
    return base_bodies;
}

} // namespace

void preload_shuttlecock_geometry() {
    (void)shuttlecockBaseBodies();
}

std::pair<Eigen::Vector3d, Eigen::Vector3d> shuttlecock_aero(double Vx, double Vy, double Vz,
                                                             double eta1_rad, double eta2_rad,
                                                             double eta3_rad, double eta4_rad,
                                                             double rho, double T_K, double s,
                                                             const ShuttlecockDebugOptions& debug_opts) {
    const auto& bodies = shuttlecockBaseBodies();

    // Apply rotation angles: choose signs so wings rotate outward about z-axis hinges
    // Order: [box, +X wing, -X wing, +Y wing, -Y wing]
    thread_local std::vector<double> angles;
    angles.assign(bodies.size(), 0.0);
    if (angles.size() >= 5) {
        angles[1] = -eta1_rad;
        angles[2] = eta2_rad;
        angles[3] = eta3_rad;
        angles[4] = -eta4_rad;
    }
    // const double pi = std::acos(-1.0);
    // const double rad2deg = 180 / pi;
    // std::cout << "WingLeft: " << angles[1]*rad2deg*-1 << '\n';
    // std::cout << "WingRight: " << angles[2]*rad2deg << '\n';
    // std::cout << "WingTop: " << angles[3]*rad2deg << '\n';
    // std::cout << "WingBottom: " << angles[4]*rad2deg*-1 << '\n';

    // Velocity and wind in I-frame
    Eigen::Vector3d V_I(Vx, Vy, Vz);
    Eigen::Vector3d wind_I = Eigen::Vector3d::Zero();
    Eigen::Quaterniond q_BI(1.0, 0.0, 0.0, 0.0); // body aligned with inertial
    Eigen::Vector3d omega_BI_B = Eigen::Vector3d::Zero();

    // Infer particle mass from s = |V| / cm, cm = sqrt(2 kB T / m)
    const double kB = 1.38064852e-23;
    double V = V_I.norm();
    if (V <= 0.0 || s <= 0.0) {
        return {Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()};
    }
    // m = 2 kB T (s / V)^2
    double m_particle = 2.0 * kB * T_K * (s / V) * (s / V);

    // Use temperature ratio method 2 by default
    int temp_ratio_method = 2;

    bool debug_vtk_enabled = debug_opts.enable_debug_vtk;
    std::string vtk_filename = debug_opts.vtk_filename.empty()
                                   ? "shuttlecock_debug.vtk"
                                   : debug_opts.vtk_filename;
    if (!debug_vtk_enabled) {
        if (const char* env = std::getenv("SHUTTLECOCK_DEBUG_VTK")) {
            debug_vtk_enabled = true;
            if (env[0] != '\0') {
                vtk_filename = env;
            }
        }
    }
    const std::string* vtk_output_ptr = debug_vtk_enabled ? &vtk_filename : nullptr;

    auto out = vleoAerodynamics(q_BI, omega_BI_B, V_I, wind_I, rho, T_K, m_particle,
                                 bodies, angles, temp_ratio_method, vtk_output_ptr);
    return {out.force, out.torque};
}

}
